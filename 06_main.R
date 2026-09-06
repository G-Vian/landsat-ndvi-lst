# =============================================================================
# MODULE 06 - MAIN ORCHESTRATOR
# Landsat Collection 2 Level-2 | NDVI & LST time series
# =============================================================================
# RUN ONLY THIS FILE. Modules 00-05 and 07-09 are sourced automatically.
#
#   Rscript R/06_main.R
#
# On an HPC cluster with a container:
#   singularity exec --bind $HOME:$HOME container.sif Rscript R/06_main.R
#
# YOU SHOULD NOT NEED TO EDIT THIS FILE to change study area, period or
# thresholds — all of that lives in 00_config.R. The few blocks here that you
# might legitimately want to change are marked ### EDIT ###.
#
# Pipeline:
#   00_config.R       -> configuration, paths, validation
#   01_scenes.R       -> scene inventory + MTL metadata
#   02_qa_mask.R      -> QA_PIXEL quality mask
#   03_calc_indices.R -> NDVI and LST
#   04_clip_aoi.R     -> reprojection and clipping to the AOI
#   05_export_stats.R -> AOI-wide statistics and export
#   07_zonal_stats.R  -> per-zone and per-pixel statistics
#   08_plots.R        -> maps and charts
#   09_anomaly_tracker.R -> anomaly report
# =============================================================================

# ---- 0. Resolve the script directory -----------------------------------------
.main_dir <- tryCatch({
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grep("--file=", args)]
  if (length(file_arg) > 0)
    dirname(normalizePath(sub("--file=", "", file_arg[1])))
  else getwd()
}, error = function(e) getwd())

if (.main_dir != getwd()) {
  setwd(.main_dir)
  message(sprintf("[06_main] Working directory set to: %s", .main_dir))
}

# ---- 1. Load modules ---------------------------------------------------------
source(file.path(.main_dir, "00_config.R"))
source(file.path(.main_dir, "01_scenes.R"))
source(file.path(.main_dir, "02_qa_mask.R"))
source(file.path(.main_dir, "03_calc_indices.R"))
source(file.path(.main_dir, "04_clip_aoi.R"))
source(file.path(.main_dir, "05_export_stats.R"))
source(file.path(.main_dir, "07_zonal_stats.R"))
source(file.path(.main_dir, "08_plots.R"))
source(file.path(.main_dir, "09_anomaly_tracker.R"))

# ---- 2. Initialise -----------------------------------------------------------
load_packages()
create_output_dirs()

# Fail fast on configuration mistakes BEFORE spending hours processing.
validate_config()

log_msg("========== PROCESSING START ==========")
log_msg(sprintf("AOI: %s | Period: %d-%d | Target CRS: %s",
                AOI_NAME, YEAR_START, YEAR_END, CRS_TARGET))

# ---- 3. Scene inventory ------------------------------------------------------
log_msg("--- STAGE 1: Scene inventory ---")
inventory <- build_inventory()

export_csv(
  inventory |> dplyr::select(-dplyr::starts_with("file_")),
  name = "scene_inventory", subdir = "scenes"
)

# ---- 4. Load the AOI (once, cached) ------------------------------------------
log_msg("--- STAGE 2: Loading area of interest ---")
aoi_boundary <- load_aoi_boundary()
aoi_zones    <- load_aoi_zones()

# ---- 5. Per-scene processing loop --------------------------------------------
log_msg("--- STAGE 3: Scene processing ---")

valid_scenes <- inventory |> dplyr::filter(complete)
log_msg(sprintf("Scenes to process: %d", nrow(valid_scenes)))

list_stats     <- vector("list", nrow(valid_scenes))
list_rasters   <- vector("list", nrow(valid_scenes))

for (i in seq_len(nrow(valid_scenes))) {

  scene <- valid_scenes[i, ]
  log_msg(sprintf("Scene %d/%d: %s", i, nrow(valid_scenes), scene$scene_id))

  # --- 5a. Footprint overlap with the AOI ---
  # If EVERY scene fails here, MIN_AOI_OVERLAP_FRAC is too high for your AOI
  # size — see 00_config.R Section 8c.
  r_test <- tryCatch(terra::rast(scene$file_lst), error = function(e) NULL)
  if (is.null(r_test)) {
    log_msg(sprintf("Could not open raster: %s", scene$scene_id), "WARN")
    next
  }
  if (!has_overlap(r_test, aoi_boundary)) {
    log_msg(sprintf("No sufficient overlap with AOI: %s", scene$scene_id), "WARN")
    rm(r_test); gc()
    next
  }
  rm(r_test); gc()

  # --- 5b. QA mask + NDVI + LST ---
  result <- process_scene(scene)
  if (is.null(result)) next

  # --- 5c. Clip to the AOI ---
  ndvi_clip <- tryCatch(clip_to_aoi(result$ndvi, aoi_boundary),
                        error = function(e) {
                          log_msg(sprintf("NDVI clip error %s: %s",
                                          scene$scene_id, e$message), "ERROR"); NULL })
  lst_clip  <- tryCatch(clip_to_aoi(result$lst, aoi_boundary),
                        error = function(e) {
                          log_msg(sprintf("LST clip error %s: %s",
                                          scene$scene_id, e$message), "ERROR"); NULL })
  if (is.null(ndvi_clip) || is.null(lst_clip)) next

  # --- 5d. AOI-specific validity check ---
  # The scene-level check in process_scene() looks at the whole scene; a scene
  # can be clear overall and completely clouded over the AOI.
  n_valid_aoi <- tryCatch(terra::global(!is.na(ndvi_clip), "sum")[[1]],
                          error = function(e) 0)
  n_total_aoi <- terra::ncell(ndvi_clip)
  pct_valid_aoi <- if (n_total_aoi > 0) 100 * n_valid_aoi / n_total_aoi else 0

  if (pct_valid_aoi < MIN_AOI_VALID_FRAC * 100) {
    log_msg(sprintf("Valid pixels over AOI below minimum (%.1f%% < %.0f%%): %s — discarding.",
                    pct_valid_aoi, MIN_AOI_VALID_FRAC * 100, scene$scene_id), "WARN")
    rm(ndvi_clip, lst_clip, result); gc()
    next
  }

  # --- 5e. NDVI/LST consistency check ---
  # A large excess of valid NDVI over valid LST pixels means the thermal band
  # is degraded for this acquisition (common in some L7 SLC-off scenes).
  n_lst_valid <- tryCatch(terra::global(!is.na(lst_clip), "sum")[[1]],
                          error = function(e) 0)
  if (n_lst_valid > 0 &&
      n_valid_aoi / max(1, n_lst_valid) > MAX_NDVI_LST_RATIO) {
    log_msg(sprintf("NDVI/LST mismatch: %d vs %d px (%.1fx): %s — discarding.",
                    n_valid_aoi, n_lst_valid,
                    n_valid_aoi / max(1, n_lst_valid), scene$scene_id), "WARN")
    rm(ndvi_clip, lst_clip, result); gc()
    next
  }

  # --- 5f. Anomaly bookkeeping (measured over the AOI only) ---
  n_at_bound   <- 0L
  n_lst_total  <- 1L   # safe denominator

  tryCatch({
    ndvi_min_a <- terra::global(ndvi_clip, "min", na.rm = TRUE)[[1]]
    ndvi_max_a <- terra::global(ndvi_clip, "max", na.rm = TRUE)[[1]]
    lst_min_a  <- terra::global(lst_clip,  "min", na.rm = TRUE)[[1]]
    lst_max_a  <- terra::global(lst_clip,  "max", na.rm = TRUE)[[1]]
    n_lst_total <- as.integer(terra::global(!is.na(lst_clip), "sum")[[1]])

    # Count pixels sitting AT the physical bounds. A 0.05 C tolerance avoids
    # flagging genuinely cold/hot pixels that merely land near the limit.
    n_below <- as.integer(dplyr::coalesce(
      terra::global(abs(lst_clip - LST_MIN) < 0.05, "sum", na.rm = TRUE)[[1]], 0))
    n_above <- as.integer(dplyr::coalesce(
      terra::global(abs(lst_clip - LST_MAX) < 0.05, "sum", na.rm = TRUE)[[1]], 0))
    n_at_bound <- n_below + n_above

    record_anomaly(
      scene_id = scene$scene_id, sensor = scene$sensor,
      data_aq = scene$acq_date, ano = scene$year, mes = scene$month,
      ndvi_raw_min = ndvi_min_a, ndvi_raw_max = ndvi_max_a, ndvi_n_desc = 0L,
      lst_raw_min = lst_min_a, lst_raw_max = lst_max_a,
      lst_n_abaixo = n_below, lst_n_acima = n_above,
      lst_n_desc = n_at_bound, lst_n_desc_dn = 0L,
      n_total_cena = as.integer(n_total_aoi)
    )

    if (n_at_bound > 0)
      log_msg(sprintf(
        "[anomaly] %s | %d px at LST_MIN (%.0fC) | %d px at LST_MAX (%.0fC) | LST [%.1f, %.1f]",
        scene$scene_id, n_below, LST_MIN, n_above, LST_MAX, lst_min_a, lst_max_a), "WARN")

  }, error = function(e)
    log_msg(sprintf("[anomaly] ERROR in %s: %s", scene$scene_id, e$message), "ERROR"))

  # --- 5g. Reject scenes whose thermal band is systematically saturated ---
  pct_at_bound <- n_at_bound / max(1L, n_lst_total)
  if (pct_at_bound > MAX_LST_AT_BOUND_FRAC) {
    log_msg(sprintf("LST has %.1f%% pixels stuck at a bound: %s — discarding.",
                    100 * pct_at_bound, scene$scene_id), "WARN")
    rm(ndvi_clip, lst_clip, result); gc()
    next
  }

  # --- 5h. Save rasters (optional; see SAVE_CLIPPED_RASTERS) ---
  if (isTRUE(SAVE_CLIPPED_RASTERS)) {
    paths <- save_rasters(ndvi_clip, lst_clip, scene$scene_id)
    list_rasters[[i]] <- list(ndvi_path = paths$ndvi,
                              lst_path  = paths$lst,
                              scene_id  = scene$scene_id)
  }

  # --- 5i. Per-scene AOI statistics ---
  scene$pct_valid  <- result$pct_valid
  list_stats[[i]]  <- extract_scene_stats(ndvi_clip, lst_clip, scene)

  rm(result, ndvi_clip, lst_clip); gc()
  log_msg(sprintf("Scene %d/%d done: %s", i, nrow(valid_scenes), scene$scene_id))
}

# ---- 6. Consolidate ----------------------------------------------------------
log_msg("--- STAGE 4: Consolidation and export ---")

df_scenes <- dplyr::bind_rows(Filter(Negate(is.null), list_stats))

if (nrow(df_scenes) == 0) {
  log_msg("NO scene processed successfully. Check the log.", "ERROR")
  stop(paste(
    "Processing finished with no results. Most common causes:",
    "  - MIN_AOI_OVERLAP_FRAC too high for a large AOI (00_config.R 8c)",
    "  - AOI shapefile does not actually intersect the downloaded scenes",
    "  - LST_MIN/LST_MAX wrong for this climate, discarding everything",
    "  - Period/tier filters excluded all scenes",
    sep = "\n"))
}

log_msg(sprintf("Scenes with usable data: %d", nrow(df_scenes)))

# ---- 7. Exports --------------------------------------------------------------
export_csv(df_scenes, "ndvi_lst_per_scene", subdir = "scenes")

df_monthly <- aggregate_monthly(df_scenes)
export_csv(df_monthly, "ndvi_lst_monthly", subdir = "monthly")

df_annual <- aggregate_annual(df_monthly)
export_csv(df_annual, "ndvi_lst_annual", subdir = "monthly")

export_txt_report(df_scenes, df_monthly)

tryCatch({
  sink(file.path(OUTPUT_DIR, "tables", "annual_summary.txt"))
  on.exit(if (sink.number() > 0) sink(), add = TRUE)
  cat(sprintf("=== ANNUAL SUMMARY - NDVI and LST - %s ===\n\n", AOI_NAME))
  print(as.data.frame(df_annual), row.names = FALSE)
  sink()
}, error = function(e) {
  if (sink.number() > 0) sink()
  log_msg(sprintf("[main] ERROR writing annual_summary.txt: %s", e$message), "ERROR")
})

# ---- 8. Zonal and pixel statistics -------------------------------------------
# Controlled by COMPUTE_ZONE_STATS / COMPUTE_PIXEL_TABLES in 00_config.R.
spatial_results <- NULL
if (isTRUE(COMPUTE_ZONE_STATS)) {
  log_msg("--- STAGE 5: Zonal statistics ---")
  spatial_results <- tryCatch(
    run_zonal_stats(df_scene_stats = df_scenes,
                    zones          = aoi_zones,
                    scene_rasters  = list_rasters),
    error = function(e) {
      log_msg(sprintf("[main] ERROR in run_zonal_stats: %s", e$message), "ERROR")
      NULL
    })
} else {
  log_msg("[main] Zonal statistics disabled (COMPUTE_ZONE_STATS = FALSE).")
}

# ---- 9. Plots ----------------------------------------------------------------
if (isTRUE(GENERATE_PLOTS) && !is.null(spatial_results)) {
  log_msg("--- STAGE 6: Maps and charts ---")
  tryCatch(
    run_plots(df_scene_stats  = df_scenes,
              df_zone_monthly = spatial_results$zone_monthly,
              df_zone_annual  = spatial_results$zone_annual,
              zones_sf        = aoi_zones,
              df_monthly      = df_monthly),
    error = function(e)
      log_msg(sprintf("[main] ERROR in run_plots: %s", e$message), "ERROR"))
} else if (!isTRUE(GENERATE_PLOTS)) {
  log_msg("[main] Plots disabled (GENERATE_PLOTS = FALSE).")
} else {
  log_msg("[main] Plots skipped: zonal statistics unavailable.", "WARN")
}

rm(list_rasters, spatial_results); gc()

# ---- 10. Anomaly report ------------------------------------------------------
log_msg("--- STAGE 7: Anomaly report ---")
tryCatch(run_anomaly_report(),
         error = function(e)
           log_msg(sprintf("[main] ERROR in anomaly report: %s", e$message), "ERROR"))

# ---- 11. Final summary -------------------------------------------------------
log_msg("========== PROCESSING COMPLETE ==========")
cat("\n=============================================================\n")
cat("                 PROCESSING FINISHED                         \n")
cat("=============================================================\n\n")
cat(sprintf("  Area of interest  : %s\n", AOI_NAME))
cat(sprintf("  Scenes processed  : %d\n", nrow(df_scenes)))
cat(sprintf("  Period covered    : %s to %s\n",
            min(df_scenes$data_aq), max(df_scenes$data_aq)))
cat(sprintf("  Mean NDVI         : %.4f\n", mean(df_scenes$ndvi_media, na.rm = TRUE)))
cat(sprintf("  Mean LST (C)      : %.2f\n", mean(df_scenes$lst_c_media, na.rm = TRUE)))
cat(sprintf("\n  Outputs in: %s\n", OUTPUT_DIR))
cat("  |-- rasters/ndvi/    per-scene NDVI GeoTIFFs\n")
cat("  |-- rasters/lst/     per-scene LST GeoTIFFs (Celsius)\n")
cat("  |-- tables/scenes/   per-scene CSV\n")
cat("  |-- tables/monthly/  monthly and annual CSV\n")
cat("  |-- tables/spatial/  per-zone and per-pixel CSV\n")
cat("  |-- tables/anomalies/ quality-control report\n")
cat("  |-- plots/           maps and charts\n")
cat("  `-- logs/            timestamped processing log\n\n")
