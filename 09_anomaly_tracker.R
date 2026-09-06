# =============================================================================
# MODULE 09 - NDVI / LST ANOMALY TRACKER
# =============================================================================
# Records, accumulates and exports a quality-control report of every pixel
# discarded or flagged as suspicious during processing.
#
# WHY THIS MATTERS WHEN ADAPTING TO A NEW AREA:
#   The anomaly report is how you VALIDATE your choice of LST_MIN / LST_MAX
#   (00_config.R Section 6). After the first run on a new region, open
#   tables/anomalies/ and check:
#
#     - Many scenes with pixels piling up AT a bound  -> your bound is cutting
#       into real data. Widen it.
#     - Nothing ever discarded, across all sensors    -> your bounds may be too
#       loose to catch genuine artefacts. Tighten them.
#     - One sensor systematically worse than others   -> expected for Landsat 7
#       after 2003 (SLC-off); consider excluding L7 or that period.
#
# Outputs under OUTPUT_DIR/tables/anomalies/:
#   anomalies_per_scene.csv / .txt   one row per scene
#   summary_by_sensor.csv   / .txt   grouped by sensor
#   summary_by_year.csv     / .txt   grouped by year
#   critical_scenes.csv     / .txt   scenes with >5% discarded pixels
#
# Integration:
#   record_anomaly()     is called from the scene loop in 06_main.R
#   run_anomaly_report() is called at the end of 06_main.R
# =============================================================================

source("00_config.R")

# =============================================================================
# 1. GLOBAL LOG - grows by one row per processed scene
# =============================================================================

.ANOMALY_LOG <- tibble::tibble(
  scene_id      = character(),
  sensor        = character(),
  data_aq       = as.Date(character()),
  ano           = integer(),
  mes           = integer(),
  ndvi_raw_min  = double(),    # NDVI range observed BEFORE clamping
  ndvi_raw_max  = double(),
  ndvi_n_desc   = integer(),   # pixels discarded by the NDVI clamp
  lst_raw_min   = double(),    # LST range observed BEFORE clamping
  lst_raw_max   = double(),
  lst_n_abaixo  = integer(),   # pixels at/below LST_MIN
  lst_n_acima   = integer(),   # pixels at/above LST_MAX
  lst_n_desc    = integer(),   # total discarded by the LST clamp
  lst_n_desc_dn = integer(),   # discarded by the sentinel-DN filter
  n_total_cena  = integer()    # denominator: total pixels considered
)

# =============================================================================
# 2. record_anomaly() - called once per scene
# =============================================================================

#' Append one scene's anomaly metrics to the global log.
#'
#' Column names are kept in the original Portuguese for compatibility with
#' existing downstream analyses; see the note in 05_export_stats.R.
record_anomaly <- function(scene_id, sensor, data_aq, ano, mes,
                           ndvi_raw_min, ndvi_raw_max, ndvi_n_desc,
                           lst_raw_min, lst_raw_max,
                           lst_n_abaixo, lst_n_acima, lst_n_desc,
                           lst_n_desc_dn, n_total_cena) {

  as_int <- function(x) if (is.null(x) || length(x) == 0 || is.na(x)) NA_integer_ else as.integer(x)
  as_dbl <- function(x) if (is.null(x) || length(x) == 0) NA_real_ else as.double(x)

  new_row <- tibble::tibble(
    scene_id      = as.character(scene_id),
    sensor        = as.character(sensor),
    data_aq       = as.Date(data_aq),
    ano           = as_int(ano),
    mes           = as_int(mes),
    ndvi_raw_min  = as_dbl(ndvi_raw_min),
    ndvi_raw_max  = as_dbl(ndvi_raw_max),
    ndvi_n_desc   = as_int(ndvi_n_desc),
    lst_raw_min   = as_dbl(lst_raw_min),
    lst_raw_max   = as_dbl(lst_raw_max),
    lst_n_abaixo  = as_int(lst_n_abaixo),
    lst_n_acima   = as_int(lst_n_acima),
    lst_n_desc    = as_int(lst_n_desc),
    lst_n_desc_dn = as_int(lst_n_desc_dn),
    n_total_cena  = as_int(n_total_cena)
  )

  .ANOMALY_LOG <<- dplyr::bind_rows(.ANOMALY_LOG, new_row)
}

# =============================================================================
# 3. ANALYSIS
# =============================================================================

#' Add derived percentages and criticality flags to the raw log.
enrich_log <- function(df) {
  df |>
    dplyr::mutate(
      pct_ndvi_desc  = dplyr::if_else(n_total_cena > 0, 100 * ndvi_n_desc   / n_total_cena, NA_real_),
      pct_lst_desc   = dplyr::if_else(n_total_cena > 0, 100 * lst_n_desc    / n_total_cena, NA_real_),
      pct_lst_abaixo = dplyr::if_else(n_total_cena > 0, 100 * lst_n_abaixo  / n_total_cena, NA_real_),
      pct_lst_acima  = dplyr::if_else(n_total_cena > 0, 100 * lst_n_acima   / n_total_cena, NA_real_),
      pct_dn_desc    = dplyr::if_else(n_total_cena > 0, 100 * lst_n_desc_dn / n_total_cena, NA_real_),

      # Flags relative to the CONFIGURED limits - these are what tell you
      # whether LST_MIN / LST_MAX are appropriate for your region.
      flag_lst_cold  = !is.na(lst_raw_min) & lst_raw_min < LST_MIN,
      flag_lst_hot   = !is.na(lst_raw_max) & lst_raw_max > LST_MAX,
      flag_ndvi_anom = (!is.na(ndvi_raw_min) & ndvi_raw_min < NDVI_MIN) |
                       (!is.na(ndvi_raw_max) & ndvi_raw_max > NDVI_MAX),

      # ### EDIT ### threshold for "critical": >5% of pixels discarded
      flag_critical  = (!is.na(pct_lst_desc) & pct_lst_desc > 5) |
                       (!is.na(pct_dn_desc)  & pct_dn_desc  > 5),

      anomaly_type = dplyr::case_when(
        lst_n_desc_dn > 0 & lst_n_desc > 0 ~ "sentinel DN + LST clamp",
        lst_n_desc_dn > 0                  ~ "sentinel DN (border/SLC-off)",
        lst_n_abaixo > 0 & lst_n_acima > 0 ~ "LST both cold and hot",
        lst_n_abaixo > 0                   ~ "LST below minimum",
        lst_n_acima  > 0                   ~ "LST above maximum",
        ndvi_n_desc  > 0                   ~ "NDVI outside [-1,1]",
        TRUE                               ~ "no anomaly"
      ),
      mes_nome = format(data_aq, "%B")
    ) |>
    dplyr::arrange(data_aq)
}

#' Per-sensor summary: reveals whether one satellite is systematically worse.
summarise_by_sensor <- function(df) {
  df |>
    dplyr::group_by(sensor) |>
    dplyr::summarise(
      n_scenes            = dplyr::n(),
      n_with_cold_lst     = sum(flag_lst_cold, na.rm = TRUE),
      n_with_hot_lst      = sum(flag_lst_hot,  na.rm = TRUE),
      n_with_sentinel_dn  = sum(lst_n_desc_dn > 0, na.rm = TRUE),
      n_critical          = sum(flag_critical, na.rm = TRUE),
      lst_raw_min_global  = min(lst_raw_min,  na.rm = TRUE),
      lst_raw_max_global  = max(lst_raw_max,  na.rm = TRUE),
      ndvi_raw_min_global = min(ndvi_raw_min, na.rm = TRUE),
      ndvi_raw_max_global = max(ndvi_raw_max, na.rm = TRUE),
      pct_lst_desc_mean   = mean(pct_lst_desc, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                ~ifelse(is.infinite(.), NA_real_, round(., 4))))
}

#' Per-year summary: temporal evolution of data quality.
summarise_by_year <- function(df) {
  df |>
    dplyr::group_by(ano) |>
    dplyr::summarise(
      n_scenes           = dplyr::n(),
      n_with_cold_lst    = sum(flag_lst_cold, na.rm = TRUE),
      n_with_hot_lst     = sum(flag_lst_hot,  na.rm = TRUE),
      n_with_sentinel_dn = sum(lst_n_desc_dn > 0, na.rm = TRUE),
      lst_raw_min_year   = min(lst_raw_min, na.rm = TRUE),
      lst_raw_max_year   = max(lst_raw_max, na.rm = TRUE),
      total_px_disc_lst  = sum(lst_n_desc, na.rm = TRUE),
      pct_lst_desc_mean  = mean(pct_lst_desc, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                ~ifelse(is.infinite(.), NA_real_, round(., 4))))
}

# =============================================================================
# 4. EXPORT
# =============================================================================

.save_csv_anom <- function(df, name, dir_out) {
  f <- file.path(dir_out, name)
  df |>
    dplyr::mutate(dplyr::across(dplyr::where(is.double), ~round(., 4))) |>
    readr::write_csv(f, na = "NA")
  log_msg(sprintf("[09_anom] CSV: %s (%d rows)", name, nrow(df)))
  invisible(f)
}

.save_txt_anom <- function(df, name, dir_out, title) {
  f <- file.path(dir_out, name)
  sink(f)
  on.exit(if (sink.number() > 0) sink(), add = TRUE)
  cat("=======================================================================\n")
  cat(sprintf("  %s | %s\n", title, AOI_NAME))
  cat(sprintf("  Generated: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("  Configured limits: NDVI [%.1f, %.1f] | LST [%.1f, %.1f] C\n",
              NDVI_MIN, NDVI_MAX, LST_MIN, LST_MAX))
  cat("=======================================================================\n\n")
  print(as.data.frame(df), row.names = FALSE)
  sink()
  log_msg(sprintf("[09_anom] TXT: %s", name))
  invisible(f)
}

# =============================================================================
# 5. ORCHESTRATOR
# =============================================================================

#' Consolidate the anomaly log and write all reports.
run_anomaly_report <- function() {

  if (nrow(.ANOMALY_LOG) == 0) {
    log_msg("[09_anom] Anomaly log is empty - nothing to report.", "WARN")
    return(invisible(NULL))
  }

  dir_out <- file.path(OUTPUT_DIR, "tables", "anomalies")
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

  enriched  <- enrich_log(.ANOMALY_LOG)
  by_sensor <- summarise_by_sensor(enriched)
  by_year   <- summarise_by_year(enriched)
  critical  <- enriched |> dplyr::filter(flag_critical)

  .save_csv_anom(enriched,  "anomalies_per_scene.csv", dir_out)
  .save_csv_anom(by_sensor, "summary_by_sensor.csv",   dir_out)
  .save_csv_anom(by_year,   "summary_by_year.csv",     dir_out)
  .save_csv_anom(critical,  "critical_scenes.csv",     dir_out)

  .save_txt_anom(
    enriched |> dplyr::select(scene_id, sensor, data_aq, lst_raw_min, lst_raw_max,
                              pct_lst_desc, anomaly_type),
    "anomalies_per_scene.txt", dir_out, "ANOMALIES PER SCENE")
  .save_txt_anom(by_sensor, "summary_by_sensor.txt", dir_out, "ANOMALY SUMMARY BY SENSOR")
  .save_txt_anom(by_year,   "summary_by_year.txt",   dir_out, "ANOMALY SUMMARY BY YEAR")
  if (nrow(critical) > 0)
    .save_txt_anom(
      critical |> dplyr::select(scene_id, sensor, data_aq, pct_lst_desc, anomaly_type),
      "critical_scenes.txt", dir_out, "CRITICAL SCENES (>5% discarded)")

  # Console guidance - the point of this module
  n_cold <- sum(enriched$flag_lst_cold, na.rm = TRUE)
  n_hot  <- sum(enriched$flag_lst_hot,  na.rm = TRUE)
  log_msg(sprintf("[09_anom] %d scenes analysed | %d with LST below LST_MIN | %d above LST_MAX | %d critical.",
                  nrow(enriched), n_cold, n_hot, nrow(critical)))

  if (n_cold > nrow(enriched) * 0.5 || n_hot > nrow(enriched) * 0.5)
    log_msg(paste("[09_anom] More than half the scenes hit an LST bound.",
                  "LST_MIN/LST_MAX in 00_config.R are probably wrong for this",
                  "climate and are discarding real data."), "WARN")

  invisible(list(per_scene = enriched, by_sensor = by_sensor,
                 by_year = by_year, critical = critical))
}
