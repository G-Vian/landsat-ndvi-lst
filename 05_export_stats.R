# =============================================================================
# MODULE 05 - AOI-WIDE STATISTICS AND EXPORT (CSV / TXT)
# =============================================================================
# Responsibility:
#   - Extract summary statistics of NDVI and LST over the WHOLE AOI, per scene
#   - Aggregate to monthly and annual scale
#   - Export CSV (machine-readable) and TXT (human-readable report)
#
# Per-ZONE statistics are the job of module 07; this module treats the AOI as
# a single unit.
#
# COLUMN NAMING NOTE
#   The statistics tables use the column names `ano`, `mes`, `data_aq`
#   (Portuguese) rather than year/month/acq_date. This is intentional: modules
#   07, 08 and 09, and several existing downstream analyses, read these names.
#   Renaming them is a breaking change - if you do it, update 07/08/09 too.
# =============================================================================

source("00_config.R")

# -----------------------------------------------------------------------------
# 1. PER-SCENE STATISTICS OVER THE AOI
# -----------------------------------------------------------------------------

#' Summary statistics of the clipped NDVI and LST rasters.
#'
#' MEMORY NOTE FOR LARGE AOIs: terra::values() pulls every pixel into RAM.
#' For a city that is trivial; for a country-sized AOI it is not. If you hit
#' memory limits, replace the body of raster_stats() with terra::global()
#' calls, which stream from disk (you lose the median, which global() does not
#' provide, but gain the ability to process arbitrarily large rasters).
#'
#' @param ndvi_r Clipped NDVI SpatRaster
#' @param lst_r  Clipped LST SpatRaster
#' @param scene  Inventory row, with $pct_valid already attached
#' @return one-row tibble
extract_scene_stats <- function(ndvi_r, lst_r, scene) {

  raster_stats <- function(r, prefix) {
    vals <- terra::values(r, na.rm = TRUE)
    if (length(vals) == 0 || all(is.na(vals))) {
      out <- tibble::tibble(a = NA_real_, b = NA_real_, c = NA_real_,
                            d = NA_real_, e = NA_real_, f = 0L)
    } else {
      out <- tibble::tibble(
        a = mean(vals,   na.rm = TRUE),
        b = median(vals, na.rm = TRUE),
        c = min(vals,    na.rm = TRUE),
        d = max(vals,    na.rm = TRUE),
        e = sd(vals,     na.rm = TRUE),
        # _n is the VALID pixel count. It doubles as the weight in every
        # weighted aggregation downstream, so never drop this column.
        f = as.integer(length(vals))
      )
    }
    names(out) <- paste0(prefix, c("_media", "_mediana", "_min", "_max", "_dp", "_n"))
    out
  }

  dplyr::bind_cols(
    tibble::tibble(
      scene_id      = scene$scene_id,
      sensor        = scene$sensor,
      data_aq       = scene$acq_date,
      ano           = scene$year,
      mes           = scene$month,
      mes_nome      = format(scene$acq_date, "%B"),
      pathrow       = scene$pathrow,
      tier          = scene$tier,
      cloud_cover   = scene$cloud_cover,
      sun_elevation = scene$sun_elevation,
      pct_valido    = scene$pct_valid
    ),
    raster_stats(ndvi_r, "ndvi"),
    raster_stats(lst_r,  "lst_c")
  )
}

# -----------------------------------------------------------------------------
# 2. TEMPORAL AGGREGATION - WEIGHTED BY VALID PIXEL COUNT
# -----------------------------------------------------------------------------
# Rationale: within a month there may be several scenes with very different
# cloud coverage. A scene contributing 40,000 clean pixels should count more
# than one contributing 400. Hence:
#
#     weighted_mean = sum(value_i * n_i) / sum(n_i)
#
# where n_i is the valid pixel count of scene i. Scenes with n_i = 0 or NA are
# excluded. This is a spatially-weighted composite, not a simple scene average.
# -----------------------------------------------------------------------------

#' Weighted mean that returns NA (rather than NaN) when all weights are zero.
weighted_mean_safe <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (sum(ok) == 0) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

# min/max that do not emit -Inf warnings on all-NA groups
safe_min <- function(x, ...) { x <- x[!is.na(x)]; if (length(x) == 0) NA_real_ else min(x) }
safe_max <- function(x, ...) { x <- x[!is.na(x)]; if (length(x) == 0) NA_real_ else max(x) }


#' Compact sensor listing for a group of scenes, e.g. "L7x2, L8x1".
#' Useful for spotting months whose composite mixes sensor generations.
sensors_compact <- function(sensors) {
  if (length(sensors) == 0 || all(is.na(sensors))) return(NA_character_)
  s_clean <- sensors[!is.na(sensors)]
  if (length(s_clean) == 0) return(NA_character_)
  s_short <- gsub("Landsat", "L", s_clean)
  tab <- table(s_short)
  tab <- tab[order(names(tab))]
  paste(sprintf("%sx%d", names(tab), as.integer(tab)), collapse = ", ")
}
# alias kept for backward compatibility with the Portuguese version
sensores_compactos <- sensors_compact


#' Merge several "L7x2, L8x1" strings into one, summing the counts.
merge_sensor_strings <- function(v) {
  lists <- v[!is.na(v)]
  if (length(lists) == 0) return(NA_character_)
  tokens <- unlist(strsplit(paste(lists, collapse = ", "), ", "))
  codes  <- sub("x.*$", "", tokens)
  counts <- suppressWarnings(as.integer(sub("^.*x", "", tokens)))
  ok <- !is.na(counts)
  if (!any(ok)) return(NA_character_)
  tab <- tapply(counts[ok], codes[ok], sum)
  tab <- tab[order(names(tab))]
  paste(sprintf("%sx%d", names(tab), as.integer(tab)), collapse = ", ")
}


#' Aggregate per-scene statistics to monthly scale.
#'
#' Produces a COMPLETE month grid over YEAR_START..YEAR_END: months with no
#' usable scene appear as explicit NA rows rather than being absent. This
#' matters for downstream time-series work, where a missing row and a zero are
#' very different things.
#'
#' CAVEAT ON THE `_dp` COLUMNS: the monthly standard deviation is a weighted
#' mean of the per-scene standard deviations, which is NOT a correct pooled SD
#' (it ignores between-scene variance). Treat it as a rough indicator of
#' within-scene spatial heterogeneity, not as a rigorous dispersion estimate.
aggregate_monthly <- function(df_scenes) {

  grid <- expand.grid(ano = YEAR_START:YEAR_END, mes = 1:12)
  grid <- tibble::as_tibble(grid)
  grid$data_ref <- as.Date(sprintf("%d-%02d-15", grid$ano, grid$mes))
  grid$mes_nome <- format(grid$data_ref, "%B")

  monthly <- df_scenes |>
    dplyr::group_by(ano, mes) |>
    dplyr::summarise(
      n_cenas         = dplyr::n(),
      sensores_usados = sensors_compact(sensor),

      # --- NDVI, weighted by valid NDVI pixels ---
      ndvi_media    = weighted_mean_safe(ndvi_media,   ndvi_n),
      ndvi_mediana  = weighted_mean_safe(ndvi_mediana, ndvi_n),
      ndvi_dp       = weighted_mean_safe(ndvi_dp,      ndvi_n),
      ndvi_min      = safe_min(ndvi_min),
      ndvi_max      = safe_max(ndvi_max),
      ndvi_n        = sum(ndvi_n, na.rm = TRUE),

      # --- LST, weighted by valid LST pixels ---
      lst_c_media   = weighted_mean_safe(lst_c_media,   lst_c_n),
      lst_c_mediana = weighted_mean_safe(lst_c_mediana, lst_c_n),
      lst_c_dp      = weighted_mean_safe(lst_c_dp,      lst_c_n),
      lst_c_min     = safe_min(lst_c_min),
      lst_c_max     = safe_max(lst_c_max),
      lst_c_n       = sum(lst_c_n, na.rm = TRUE),

      # --- metadata ---
      cloud_cover_med = weighted_mean_safe(cloud_cover,   ndvi_n),
      sun_elev_med    = weighted_mean_safe(sun_elevation, ndvi_n),

      .groups = "drop"
    ) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                ~ifelse(is.infinite(.), NA_real_, .)))

  result <- grid |>
    dplyr::left_join(monthly, by = c("ano", "mes")) |>
    dplyr::arrange(data_ref)

  log_msg(sprintf("[05_export] Monthly: %d months with data, %d empty (NA).",
                  sum(!is.na(result$ndvi_media)), sum(is.na(result$ndvi_media))))
  result
}


#' Aggregate the monthly table to annual scale, weighting each month by its
#' valid pixel count (so a month represented by a single clear scene does not
#' count the same as a month with full coverage).
aggregate_annual <- function(df_monthly) {
  df_monthly |>
    dplyr::group_by(ano) |>
    dplyr::summarise(
      n_meses_com_dados = sum(!is.na(ndvi_media)),
      sensores_usados   = merge_sensor_strings(sensores_usados),
      ndvi_media_anual  = weighted_mean_safe(ndvi_media, ndvi_n),
      ndvi_min_anual    = safe_min(ndvi_min),
      ndvi_max_anual    = safe_max(ndvi_max),
      ndvi_n_total      = sum(ndvi_n, na.rm = TRUE),
      lst_media_anual   = weighted_mean_safe(lst_c_media, lst_c_n),
      lst_min_anual     = safe_min(lst_c_min),
      lst_max_anual     = safe_max(lst_c_max),
      lst_c_n_total     = sum(lst_c_n, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                ~ifelse(is.infinite(.), NA_real_, .)))
}

# -----------------------------------------------------------------------------
# 3. CSV EXPORT
# -----------------------------------------------------------------------------

#' Write a table to OUTPUT_DIR/tables/<subdir>/<name>.csv
#'
#' Uses comma-separated, dot-decimal (readr::write_csv). If your locale expects
#' semicolon/comma-decimal for Excel, switch to readr::write_csv2 here - but
#' note that any downstream script reading these files must match.
export_csv <- function(df, name, subdir = "scenes") {

  dir_out <- file.path(OUTPUT_DIR, "tables", subdir)
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

  out_file <- file.path(dir_out, paste0(name, ".csv"))

  df_out <- df |>
    dplyr::mutate(dplyr::across(dplyr::where(is.double), ~round(., 6)))

  readr::write_csv(df_out, out_file, na = "NA")
  log_msg(sprintf("[05_export] CSV saved: %s (%d rows)",
                  basename(out_file), nrow(df_out)))
  invisible(out_file)
}

# -----------------------------------------------------------------------------
# 4. TXT REPORT
# -----------------------------------------------------------------------------

#' Human-readable summary report of scenes and monthly series.
export_txt_report <- function(df_scenes, df_monthly) {

  dir_out  <- file.path(OUTPUT_DIR, "tables")
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(dir_out, "full_report.txt")

  sink(out_file)
  on.exit(if (sink.number() > 0) sink(), add = TRUE)

  cat("=======================================================================\n")
  cat(sprintf("  NDVI AND LST REPORT - %s - Landsat Collection 2 L2SP\n", AOI_NAME))
  cat(sprintf("  Generated: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("  Period: %d to %d\n", YEAR_START, YEAR_END))
  cat("  Monthly aggregation: WEIGHTED MEAN by valid pixel count\n")
  cat("=======================================================================\n\n")

  cat("--- OVERALL SUMMARY ---\n")
  cat(sprintf("  Scenes processed  : %d\n", nrow(df_scenes)))
  cat(sprintf("  Sensors used      : %s\n",
              paste(unique(df_scenes$sensor), collapse = ", ")))
  cat(sprintf("  Mean NDVI         : %.4f\n", mean(df_scenes$ndvi_media, na.rm = TRUE)))
  cat(sprintf("  Mean LST (C)      : %.2f\n", mean(df_scenes$lst_c_media, na.rm = TRUE)))
  cat("\n")

  fmt_val <- function(x, fmt) ifelse(is.na(x), "      NA", sprintf(fmt, x))

  cat("--- PER SCENE ---\n")
  cat(sprintf("%-45s %-11s %9s %9s %8s %8s\n",
              "Scene ID", "Date", "NDVI_avg", "LST_avg", "%Valid", "Cloud%"))
  cat(strrep("-", 95), "\n")
  for (i in seq_len(nrow(df_scenes))) {
    r <- df_scenes[i, ]
    cat(sprintf("%-45s %-11s %9s %9s %8s %8s\n",
                r$scene_id, as.character(r$data_aq),
                fmt_val(r$ndvi_media,  "%9.4f"),
                fmt_val(r$lst_c_media, "%9.2f"),
                fmt_val(r$pct_valido,  "%8.1f"),
                fmt_val(r$cloud_cover, "%8.1f")))
  }
  cat("\n")

  cat("--- MONTHLY SERIES ---\n")
  cat(sprintf("%-6s %-4s %-11s %-9s %-9s %-9s %-9s %-9s %-9s\n",
              "Year", "Mon", "Month",
              "NDVI_avg", "NDVI_min", "NDVI_max",
              "LST_avg", "LST_min", "LST_max"))
  cat(strrep("-", 95), "\n")
  for (i in seq_len(nrow(df_monthly))) {
    r <- df_monthly[i, ]
    f <- function(x, fmt = "%9.4f") ifelse(is.na(x), "       NA", sprintf(fmt, x))
    cat(sprintf("%-6d %-4d %-11s %-9s %-9s %-9s %-9s %-9s %-9s\n",
                r$ano, r$mes, r$mes_nome,
                f(r$ndvi_media), f(r$ndvi_min), f(r$ndvi_max),
                f(r$lst_c_media, "%9.2f"), f(r$lst_c_min, "%9.2f"), f(r$lst_c_max, "%9.2f")))
  }
  sink()

  log_msg(sprintf("[05_export] TXT report saved: %s", out_file))
  invisible(out_file)
}
