# =============================================================================
# MODULE 07 - ZONAL AND PER-PIXEL STATISTICS
# =============================================================================
# (Formerly 07_stats_espaciais.R, specialised for Santos neighbourhoods.
#  Now generic: a "zone" is any polygon in your AOI file.)
#
# Responsibility:
#   - Aggregate NDVI and LST at two spatial scales:
#       (1) PER ZONE  -> one row per polygon of SHP_FILE
#       (2) PER PIXEL -> one row per 30 x 30 m cell   [optional, see below]
#   - At two temporal scales: monthly and annual
#   - Export CSV + human-readable TXT for both
#
# Outputs under OUTPUT_DIR/tables/spatial/:
#   zones/   ndvi_lst_per_zone_monthly.csv / .txt
#            ndvi_lst_per_zone_annual.csv  / .txt
#   pixels/  ndvi_per_pixel_monthly.csv, lst_per_pixel_monthly.csv
#            ndvi_per_pixel_annual.csv,  lst_per_pixel_annual.csv
#
# ######################### SCALING WARNING #########################
# The PER-PIXEL tables have one row per pixel per period. For a small city
# (~40 km²) that is ~44,000 rows/month — fine. For a large state it is
# hundreds of millions of rows and will exhaust memory and disk.
# They are therefore controlled by COMPUTE_PIXEL_TABLES in 00_config.R and
# default to FALSE. Zonal statistics scale with the number of polygons, not
# the area, and are safe to leave on.
# ##################################################################
#
# Weighted-average convention (same as module 05): within a period, scenes
# with more valid pixels carry proportionally more weight.
# =============================================================================

source("00_config.R")
source("04_clip_aoi.R")
source("05_export_stats.R")   # reuses weighted_mean_safe(), safe_min/max()

# =============================================================================
# SECTION A - PER-ZONE STATISTICS
# =============================================================================

# -----------------------------------------------------------------------------
# A1. ZONAL EXTRACTION FOR ONE SCENE
# -----------------------------------------------------------------------------

#' Extract NDVI and LST statistics for every zone, for one scene.
#'
#' NOTE ON ZONES WITH FEW PIXELS: after QA masking (and especially with the
#' water bit enabled) a small or water-dominated zone can end up with a handful
#' of valid pixels, or none. Those rows get *_n = 0 and NA statistics. ALWAYS
#' check the *_n columns before interpreting a zone mean — a mean over 3 pixels
#' is not comparable to a mean over 30,000.
#'
#' @param ndvi_r  Clipped NDVI SpatRaster
#' @param lst_r   Clipped LST SpatRaster
#' @param zones   sf object from load_aoi_zones() (has a `zone_name` column)
#' @param scene   One row of the per-scene statistics table
#' @return tibble: one row per zone
extract_zone_stats <- function(ndvi_r, lst_r, zones, scene) {

  sid <- scene$scene_id
  zv  <- terra::vect(zones)

  zonal_stats <- function(r, prefix) {

    # terra::extract() with fun = NULL returns every pixel of every polygon,
    # tagged with the polygon ID. That is what lets us compute the median,
    # which terra::zonal() cannot do.
    ext_df <- terra::extract(r, zv, fun = NULL, na.rm = TRUE, ID = TRUE)
    colnames(ext_df)[2] <- "value"
    ext_df <- ext_df[!is.na(ext_df$value), ]

    if (nrow(ext_df) == 0) {
      return(tibble::tibble(
        zone_id = seq_len(nrow(zones)),
        !!paste0(prefix, "_media")   := NA_real_,
        !!paste0(prefix, "_mediana") := NA_real_,
        !!paste0(prefix, "_min")     := NA_real_,
        !!paste0(prefix, "_max")     := NA_real_,
        !!paste0(prefix, "_dp")      := NA_real_,
        !!paste0(prefix, "_n")       := 0L
      ))
    }

    ext_df |>
      dplyr::group_by(ID) |>
      dplyr::summarise(
        !!paste0(prefix, "_media")   := mean(value,   na.rm = TRUE),
        !!paste0(prefix, "_mediana") := median(value, na.rm = TRUE),
        !!paste0(prefix, "_min")     := min(value,    na.rm = TRUE),
        !!paste0(prefix, "_max")     := max(value,    na.rm = TRUE),
        !!paste0(prefix, "_dp")      := sd(value,     na.rm = TRUE),
        !!paste0(prefix, "_n")       := dplyr::n(),
        .groups = "drop"
      ) |>
      dplyr::rename(zone_id = ID)
  }

  stats_ndvi <- zonal_stats(ndvi_r, "ndvi")
  stats_lst  <- zonal_stats(lst_r,  "lst_c")

  id_zones <- tibble::tibble(
    zone_id   = seq_len(nrow(zones)),
    zone_name = sf::st_drop_geometry(zones)[["zone_name"]]
  )

  id_zones |>
    dplyr::left_join(stats_ndvi, by = "zone_id") |>
    dplyr::left_join(stats_lst,  by = "zone_id") |>
    dplyr::mutate(
      scene_id   = sid,
      sensor     = scene$sensor,
      data_aq    = scene$data_aq,
      ano        = scene$ano,
      mes        = scene$mes,
      mes_nome   = format(scene$data_aq, "%B"),
      pct_valido = scene$pct_valido,
      .before    = 1
    )
}

# -----------------------------------------------------------------------------
# A2. MONTHLY AGGREGATION PER ZONE
# -----------------------------------------------------------------------------

#' Aggregate per-zone scene statistics to monthly scale.
#'
#' Emits a COMPLETE zone x year x month grid: zone-months with no usable scene
#' appear as explicit NA rows. For a long series this makes the gaps visible
#' instead of silently shortening the series.
aggregate_zone_monthly <- function(df_zone) {

  unique_zones <- df_zone |> dplyr::distinct(zone_id, zone_name)

  grid <- tidyr::expand_grid(
    zone_id = unique_zones$zone_id,
    ano     = YEAR_START:YEAR_END,
    mes     = 1:12
  ) |>
    dplyr::left_join(unique_zones, by = "zone_id") |>
    dplyr::mutate(
      data_ref = as.Date(sprintf("%d-%02d-15", ano, mes)),
      mes_nome = format(data_ref, "%B")
    )

  monthly <- df_zone |>
    dplyr::group_by(zone_id, zone_name, ano, mes) |>
    dplyr::summarise(
      n_cenas         = dplyr::n(),
      sensores_usados = sensors_compact(sensor),

      ndvi_media   = weighted_mean_safe(ndvi_media,   ndvi_n),
      ndvi_mediana = weighted_mean_safe(ndvi_mediana, ndvi_n),
      ndvi_dp      = weighted_mean_safe(ndvi_dp,      ndvi_n),
      ndvi_min     = safe_min(ndvi_min),
      ndvi_max     = safe_max(ndvi_max),
      ndvi_n       = sum(ndvi_n, na.rm = TRUE),

      lst_c_media   = weighted_mean_safe(lst_c_media,   lst_c_n),
      lst_c_mediana = weighted_mean_safe(lst_c_mediana, lst_c_n),
      lst_c_dp      = weighted_mean_safe(lst_c_dp,      lst_c_n),
      lst_c_min     = safe_min(lst_c_min),
      lst_c_max     = safe_max(lst_c_max),
      lst_c_n       = sum(lst_c_n, na.rm = TRUE),

      .groups = "drop"
    ) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                ~ifelse(is.infinite(.), NA_real_, .)))

  grid |>
    dplyr::left_join(monthly, by = c("zone_id", "zone_name", "ano", "mes")) |>
    dplyr::arrange(zone_id, data_ref)
}

# -----------------------------------------------------------------------------
# A3. ANNUAL AGGREGATION PER ZONE
# -----------------------------------------------------------------------------

aggregate_zone_annual <- function(df_monthly) {
  df_monthly |>
    dplyr::group_by(zone_id, zone_name, ano) |>
    dplyr::summarise(
      n_meses_com_dados = sum(!is.na(ndvi_media)),
      sensores_usados = {
        lists <- sensores_usados[!is.na(sensores_usados)]
        if (length(lists) == 0) NA_character_ else {
          tokens <- unlist(strsplit(paste(lists, collapse = ", "), ", "))
          codes  <- sub("x.*$", "", tokens)
          counts <- as.integer(sub("^.*x", "", tokens))
          tab <- tapply(counts, codes, sum)
          tab <- tab[order(names(tab))]
          paste(sprintf("%sx%d", names(tab), as.integer(tab)), collapse = ", ")
        }
      },
      ndvi_media_anual = weighted_mean_safe(ndvi_media, ndvi_n),
      ndvi_min_anual   = safe_min(ndvi_min),
      ndvi_max_anual   = safe_max(ndvi_max),
      ndvi_n_total     = sum(ndvi_n, na.rm = TRUE),
      lst_media_anual  = weighted_mean_safe(lst_c_media, lst_c_n),
      lst_min_anual    = safe_min(lst_c_min),
      lst_max_anual    = safe_max(lst_c_max),
      lst_c_n_total    = sum(lst_c_n, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                ~ifelse(is.infinite(.), NA_real_, .))) |>
    dplyr::arrange(zone_id, ano)
}

# -----------------------------------------------------------------------------
# A4. EXPORT - ZONES
# -----------------------------------------------------------------------------

export_zone_tables <- function(df_monthly, df_annual) {

  dir_out <- file.path(OUTPUT_DIR, "tables", "spatial", "zones")
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

  readr::write_csv(
    df_monthly |> dplyr::mutate(dplyr::across(dplyr::where(is.double), ~round(., 6))),
    file.path(dir_out, "ndvi_lst_per_zone_monthly.csv"), na = "NA")
  readr::write_csv(
    df_annual |> dplyr::mutate(dplyr::across(dplyr::where(is.double), ~round(., 6))),
    file.path(dir_out, "ndvi_lst_per_zone_annual.csv"), na = "NA")

  # --- readable monthly report ---
  sink(file.path(dir_out, "ndvi_lst_per_zone_monthly.txt"))
  on.exit(if (sink.number() > 0) sink(), add = TRUE)
  cat("=======================================================================\n")
  cat(sprintf("  NDVI AND LST BY %s - MONTHLY | %s | Landsat Collection 2\n",
              toupper(ZONE_LABEL), AOI_NAME))
  cat(sprintf("  Generated: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat("  Weighted mean by valid pixel count. NA = no usable scene.\n")
  cat("  Always check *_n (pixel count) before interpreting a mean.\n")
  cat("=======================================================================\n\n")
  cat(sprintf("%-28s %-6s %-4s %10s %10s %10s\n",
              toupper(ZONE_LABEL), "Year", "Mon", "NDVI_avg", "LST_avg", "N_px"))
  cat(strrep("-", 74), "\n")
  f <- function(x, fmt) ifelse(is.na(x), "        NA", sprintf(fmt, x))
  for (i in seq_len(nrow(df_monthly))) {
    r <- df_monthly[i, ]
    cat(sprintf("%-28s %-6d %-4d %10s %10s %10s\n",
                substr(r$zone_name, 1, 28), r$ano, r$mes,
                f(r$ndvi_media, "%10.4f"), f(r$lst_c_media, "%10.2f"),
                f(r$ndvi_n, "%10.0f")))
  }
  sink()

  # --- readable annual report ---
  sink(file.path(dir_out, "ndvi_lst_per_zone_annual.txt"))
  cat("=======================================================================\n")
  cat(sprintf("  NDVI AND LST BY %s - ANNUAL | %s | Landsat Collection 2\n",
              toupper(ZONE_LABEL), AOI_NAME))
  cat("=======================================================================\n\n")
  cat(sprintf("%-28s %-6s %10s %10s %8s\n",
              toupper(ZONE_LABEL), "Year", "NDVI_avg", "LST_avg", "Months"))
  cat(strrep("-", 66), "\n")
  for (i in seq_len(nrow(df_annual))) {
    r <- df_annual[i, ]
    cat(sprintf("%-28s %-6d %10s %10s %8d\n",
                substr(r$zone_name, 1, 28), r$ano,
                f(r$ndvi_media_anual, "%10.4f"), f(r$lst_media_anual, "%10.2f"),
                r$n_meses_com_dados))
  }
  sink()

  log_msg(sprintf("[07_zonal] Zone tables exported to %s", dir_out))
}

# =============================================================================
# SECTION B - PER-PIXEL STATISTICS  (optional; see COMPUTE_PIXEL_TABLES)
# =============================================================================

# -----------------------------------------------------------------------------
# B1. WEIGHTED PIXELWISE MEAN
# -----------------------------------------------------------------------------

#' Weighted mean across a list of aligned rasters, computed pixel by pixel.
#'
#' IMPORTANT CORRECTION vs. earlier versions of this pipeline:
#' the denominator must accumulate ONLY the weights of scenes that actually
#' have a valid value AT THAT PIXEL. Summing all weights unconditionally makes
#' every pixel that is valid in only some scenes come out systematically too
#' low (it divides a partial sum by the full weight). The denominator is
#' therefore itself a raster, not a scalar.
#'
#' @param rasters list of SpatRasters
#' @param weights numeric vector, same length (valid pixel count per scene)
#' @return SpatRaster of the weighted mean, NA where no scene had data
pixel_weighted_mean <- function(rasters, weights) {

  n <- length(rasters)
  if (n == 0) return(NULL)
  if (n == 1) return(rasters[[1]])

  template <- rasters[[1]]
  for (k in seq_along(rasters)) {
    if (!terra::compareGeom(rasters[[k]], template, stopOnError = FALSE))
      rasters[[k]] <- terra::resample(rasters[[k]], template, method = "bilinear")
  }

  num <- template * 0        # weighted sum of values
  den <- template * 0        # weighted sum of CONTRIBUTING weights (per pixel)

  for (k in seq_along(rasters)) {
    r_k <- rasters[[k]]
    w_k <- weights[k]
    valid <- !is.na(r_k)
    num <- num + terra::ifel(valid, r_k * w_k, 0)
    den <- den + terra::ifel(valid, w_k,       0)
  }

  # Pixels with no contributing scene get NA rather than a division by zero
  terra::ifel(den > 0, num / den, NA)
}

# -----------------------------------------------------------------------------
# B2. COLLECT SAVED RASTERS FOR A GROUP OF SCENES
# -----------------------------------------------------------------------------

#' Load the clipped GeoTIFFs for a set of scenes, with their weights.
#' Requires SAVE_CLIPPED_RASTERS = TRUE in a previous or the current run.
collect_rasters <- function(df_scenes, dir_ndvi, dir_lst) {

  r_ndvi <- list(); r_lst <- list(); weights <- numeric(0)

  for (i in seq_len(nrow(df_scenes))) {
    sid <- df_scenes$scene_id[i]
    f_n <- file.path(dir_ndvi, paste0(sid, "_NDVI.tif"))
    f_l <- file.path(dir_lst,  paste0(sid, "_LST_Celsius.tif"))

    if (!file.exists(f_n) || !file.exists(f_l)) {
      log_msg(sprintf("[07_zonal] Raster missing, skipping: %s", sid), "WARN")
      next
    }

    w <- df_scenes$ndvi_n[i]
    if (is.null(w) || is.na(w)) w <- 1

    r_ndvi  <- c(r_ndvi, list(terra::rast(f_n)))
    r_lst   <- c(r_lst,  list(terra::rast(f_l)))
    weights <- c(weights, w)
  }

  list(ndvi = r_ndvi, lst = r_lst, weights = weights)
}

# -----------------------------------------------------------------------------
# B3. PIXEL TABLES
# -----------------------------------------------------------------------------

#' Convert a raster to a long data frame of (x, y, value), dropping NA cells.
raster_to_pixel_df <- function(r, col_val = "valor") {
  df <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  if (ncol(df) < 3) return(tibble::tibble())
  colnames(df)[3] <- col_val
  tibble::as_tibble(df)
}

#' Build per-pixel tables for each period.
#'
#' ### SCALING: see the warning at the top of this file. ###
build_pixel_tables <- function(df_scene_stats, scale = "monthly",
                               dir_ndvi, dir_lst) {

  if (scale == "monthly") {
    groups <- df_scene_stats |>
      dplyr::group_by(ano, mes) |> dplyr::group_split() |>
      lapply(function(g) list(label = sprintf("%04d-%02d", g$ano[1], g$mes[1]),
                              ano = g$ano[1], mes = g$mes[1], scenes = g))
  } else {
    groups <- df_scene_stats |>
      dplyr::group_by(ano) |> dplyr::group_split() |>
      lapply(function(g) list(label = sprintf("%04d", g$ano[1]),
                              ano = g$ano[1], mes = NA_integer_, scenes = g))
  }

  out_ndvi <- list(); out_lst <- list()

  for (g in groups) {
    rs <- collect_rasters(g$scenes, dir_ndvi, dir_lst)
    if (length(rs$ndvi) == 0) next

    m_ndvi <- pixel_weighted_mean(rs$ndvi, rs$weights)
    m_lst  <- pixel_weighted_mean(rs$lst,  rs$weights)

    if (!is.null(m_ndvi)) {
      d <- raster_to_pixel_df(m_ndvi, "ndvi")
      if (nrow(d) > 0)
        out_ndvi[[length(out_ndvi) + 1]] <-
          dplyr::mutate(d, ano = g$ano, mes = g$mes, periodo = g$label, .before = 1)
    }
    if (!is.null(m_lst)) {
      d <- raster_to_pixel_df(m_lst, "lst_c")
      if (nrow(d) > 0)
        out_lst[[length(out_lst) + 1]] <-
          dplyr::mutate(d, ano = g$ano, mes = g$mes, periodo = g$label, .before = 1)
    }
    rm(rs, m_ndvi, m_lst); gc()
  }

  list(ndvi = dplyr::bind_rows(out_ndvi), lst = dplyr::bind_rows(out_lst))
}

export_pixel_tables <- function(tbl_monthly, tbl_annual) {

  dir_out <- file.path(OUTPUT_DIR, "tables", "spatial", "pixels")
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

  wr <- function(df, name) {
    if (is.null(df) || nrow(df) == 0) {
      log_msg(sprintf("[07_zonal] Pixel table empty, not written: %s", name), "WARN")
      return(invisible(NULL))
    }
    p <- file.path(dir_out, name)
    readr::write_csv(df |> dplyr::mutate(dplyr::across(dplyr::where(is.double), ~round(., 6))),
                     p, na = "NA")
    log_msg(sprintf("[07_zonal] Pixel table saved: %s (%d rows)", name, nrow(df)))
  }

  wr(tbl_monthly$ndvi, "ndvi_per_pixel_monthly.csv")
  wr(tbl_monthly$lst,  "lst_per_pixel_monthly.csv")
  wr(tbl_annual$ndvi,  "ndvi_per_pixel_annual.csv")
  wr(tbl_annual$lst,   "lst_per_pixel_annual.csv")
}

# =============================================================================
# ORCHESTRATOR
# =============================================================================

#' Run all zonal (and optionally per-pixel) statistics.
#'
#' @param df_scene_stats per-scene statistics table (module 05)
#' @param zones          sf zones from load_aoi_zones()
#' @param scene_rasters  list of list(ndvi_path, lst_path, scene_id)
#' @return list(zone_monthly, zone_annual, pixel_monthly, pixel_annual)
run_zonal_stats <- function(df_scene_stats, zones, scene_rasters) {

  log_msg("[07_zonal] === START - Spatial statistics ===")

  # --- SECTION A: per zone ---
  log_msg(sprintf("[07_zonal] A - Extracting per-%s statistics scene by scene...",
                  ZONE_LABEL))

  # Look up scene metadata BY scene_id, never by position: scene_rasters is
  # indexed over all candidate scenes (with NULL holes for skipped ones) while
  # df_scene_stats only holds the successful ones, so the indices differ.
  stats_lookup <- split(df_scene_stats, seq_len(nrow(df_scene_stats)))
  names(stats_lookup) <- df_scene_stats$scene_id

  list_zone <- vector("list", length(scene_rasters))

  for (i in seq_along(scene_rasters)) {
    item <- scene_rasters[[i]]
    if (is.null(item)) next
    scene <- stats_lookup[[item$scene_id]]
    if (is.null(scene)) next

    tryCatch({
      list_zone[[i]] <- extract_zone_stats(
        ndvi_r = terra::rast(item$ndvi_path),
        lst_r  = terra::rast(item$lst_path),
        zones  = zones,
        scene  = scene)
    }, error = function(e)
      log_msg(sprintf("[07_zonal] ERROR on scene %s: %s",
                      item$scene_id, conditionMessage(e)), "ERROR"))
  }

  df_zone_scenes <- dplyr::bind_rows(list_zone[!sapply(list_zone, is.null)])

  if (nrow(df_zone_scenes) == 0) {
    log_msg("[07_zonal] No zonal statistics could be extracted.", "WARN")
    return(invisible(NULL))
  }

  log_msg(sprintf("[07_zonal] Zonal rows compiled: %d (%d zones x scenes).",
                  nrow(df_zone_scenes), length(unique(df_zone_scenes$zone_id))))

  df_zone_monthly <- aggregate_zone_monthly(df_zone_scenes)
  df_zone_annual  <- aggregate_zone_annual(df_zone_monthly)
  export_zone_tables(df_zone_monthly, df_zone_annual)

  # --- SECTION B: per pixel (optional) ---
  tbl_monthly <- list(ndvi = tibble::tibble(), lst = tibble::tibble())
  tbl_annual  <- list(ndvi = tibble::tibble(), lst = tibble::tibble())

  if (isTRUE(COMPUTE_PIXEL_TABLES)) {
    log_msg("[07_zonal] B - Building per-pixel tables...")
    dir_ndvi <- file.path(OUTPUT_DIR, "rasters", "ndvi")
    dir_lst  <- file.path(OUTPUT_DIR, "rasters", "lst")
    tbl_monthly <- build_pixel_tables(df_scene_stats, "monthly", dir_ndvi, dir_lst)
    tbl_annual  <- build_pixel_tables(df_scene_stats, "annual",  dir_ndvi, dir_lst)
    export_pixel_tables(tbl_monthly, tbl_annual)
  } else {
    log_msg("[07_zonal] B - Per-pixel tables disabled (COMPUTE_PIXEL_TABLES = FALSE).")
  }

  log_msg("[07_zonal] === END - Spatial statistics ===")

  invisible(list(zone_monthly  = df_zone_monthly,
                 zone_annual   = df_zone_annual,
                 pixel_monthly = tbl_monthly,
                 pixel_annual  = tbl_annual))
}

# --- backward-compatibility alias ---
run_stats_espaciais <- run_zonal_stats
