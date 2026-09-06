# =============================================================================
# VALIDATION - 01 - NDVI / LST COVERAGE CALENDAR (2010-2025)
# =============================================================================
# QUESTION ANSWERED
#   "Which months of the time series actually carry NDVI and LST data, and
#    where are the gaps?"
#
# WHY THIS MATTERS
#   Landsat revisits a given path every 16 days, and a scene is only usable if
#   it is not obscured by cloud. In a humid subtropical coastal city such as
#   many tropical regions the rainy season routinely removes entire months from
#   record. Any statement about an annual mean therefore has to be read
#   together with how many months contributed to it: an "annual mean" built
#   from three scattered summer scenes is not comparable to one built from
#   eleven months spread across the seasons. This script makes that structure
#   explicit before any inference is drawn from the series.
#
# WHAT THE SCRIPT DOES
#   1. Reads the monthly tables written by the pipeline (city-wide and
#      per-zone), locating them by column content rather than by file name.
#   2. Builds a 16-year x 12-month calendar and classifies each cell as
#      carrying both variables, only one, or nothing.
#   3. At zone level, computes the fraction of zones with data in each month,
#      which exposes partial months where a scene exists but covers only part
#      of the city.
#   4. Flags any month declared low quality in 00_common.R separately.
#   5. Writes CSV tables, a plain-text report and two heat-map figures.
#
# HOW TO RUN
#   singularity exec --bind $HOME:$HOME \
#     /path/to/container.sif \
#     Rscript 01_coverage_calendar.R
# =============================================================================

# --- locate 00_common.R -------------------------------------------------------
# When launched as `Rscript /full/path/01_....R` the working directory is
# wherever the shell happened to be, not where the script lives, so a bare
# source("00_common.R") would fail. This resolves the script's own directory
# from the command line and falls back to the working directory when running
# interactively.
.this_dir <- tryCatch({
  a <- commandArgs(trailingOnly = FALSE)
  f <- a[grep("--file=", a)]
  if (length(f) > 0) dirname(normalizePath(sub("--file=", "", f[1])))
  else getwd()
}, error = function(e) getwd())

.common <- file.path(.this_dir, "00_common.R")
if (!file.exists(.common)) .common <- file.path(getwd(), "00_common.R")
if (!file.exists(.common))
  stop("00_common.R not found next to this script or in the working directory.\n",
       "  Looked in: ", .this_dir, "\n  and       : ", getwd())
source(.common)

TAG     <- "[01_coverage]"
DIR_OUT <- file.path(VALIDATION_ROOT, "01_coverage_calendar")
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

message(sprintf("\n%s Start: %s", TAG, Sys.time()))
message(sprintf("%s Output: %s", TAG, DIR_OUT))


# -----------------------------------------------------------------------------
# 1. LOAD DATA
# -----------------------------------------------------------------------------
pipeline_dir <- resolve_pipeline_dir()
tables_dir   <- resolve_tables_dir(pipeline_dir)
message(sprintf("%s Pipeline results: %s", TAG, pipeline_dir))

city <- load_city_monthly(tables_dir)

# Resolve the study period. YEAR_START / YEAR_END may be left NULL in
# 00_common.R, in which case the span is taken from the data itself, so the
# report always describes the run that actually happened.
resolve_period(city)

# The per-zone table is optional: without it the script still produces the
# city-level calendar, just not the spatial-coverage panel.
zones <- tryCatch(load_zone_monthly(tables_dir), error = function(e) {
  message(TAG, " Per-zone table unavailable, skipping the zone panel.")
  NULL
})


# -----------------------------------------------------------------------------
# 2. CITY-LEVEL CALENDAR
# -----------------------------------------------------------------------------
# expand.grid guarantees that months with NO scene at all appear as rows. A
# left join onto the pipeline output alone would silently drop them, and the
# gaps - which are the whole point of this figure - would become invisible.

keep_cols <- intersect(c("year", "month", "ndvi_mean", "lst_mean",
                         "n_scenes", "sensors_used"), names(city))

calendar <- expand.grid(year  = YEAR_START:YEAR_END,
                        month = 1:12) |>
  dplyr::as_tibble() |>
  dplyr::left_join(dplyr::select(city, dplyr::all_of(keep_cols)),
                   by = c("year", "month")) |>
  dplyr::mutate(
    has_ndvi = !is.na(ndvi_mean),
    has_lst  = !is.na(lst_mean),
    status   = dplyr::case_when(
      has_ndvi &  has_lst ~ "NDVI + LST",
      has_ndvi & !has_lst ~ "NDVI only",
      !has_ndvi & has_lst ~ "LST only",
      TRUE                ~ "No data"
    )
  )

# Overwrite the status of flagged months, but only where data actually exists:
# a month that is empty anyway should stay labelled "No data" rather than being
# promoted to a quality warning.
# The join is skipped entirely when no months are flagged, which is the default
# state of a fresh checkout.
if (nrow(LOW_QUALITY_MONTHS) > 0) {
  calendar <- calendar |>
    dplyr::left_join(
      dplyr::transmute(LOW_QUALITY_MONTHS, year, month, flagged = TRUE),
      by = c("year", "month")
    ) |>
    dplyr::mutate(
      status = ifelse(!is.na(flagged) & status != "No data",
                      "Flagged: low quality", status)
    )
}


# -----------------------------------------------------------------------------
# 3. ZONE-LEVEL COVERAGE
# -----------------------------------------------------------------------------
# A month can have a scene and still leave zones empty: cloud, or the SLC-off
# stripes of Landsat 7, can fall over part of the city. The percentage below is
# the share of zones with at least one valid pixel, which is a stricter and
# more honest measure of coverage than "a scene exists".

zone_cover <- NULL
if (!is.null(zones)) {
  n_zones_total <- dplyr::n_distinct(zones$zone_name)

  zone_cover <- zones |>
    dplyr::group_by(year, month) |>
    dplyr::summarise(
      n_zones_lst  = sum(!is.na(lst_mean)),
      n_zones_ndvi = sum(!is.na(ndvi_mean)),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      n_zones_total = n_zones_total,
      pct_lst  = round(100 * n_zones_lst  / n_zones_total, 1),
      pct_ndvi = round(100 * n_zones_ndvi / n_zones_total, 1)
    )

  # Fill in months with no scene at all as 0 per cent rather than NA, so the
  # heat map shows an explicit empty cell instead of a hole.
  zone_cover <- expand.grid(year = YEAR_START:YEAR_END, month = 1:12) |>
    dplyr::as_tibble() |>
    dplyr::left_join(zone_cover, by = c("year", "month")) |>
    dplyr::mutate(
      dplyr::across(dplyr::starts_with("pct"),  ~ifelse(is.na(.), 0,  .)),
      dplyr::across(dplyr::starts_with("n_zones"), ~ifelse(is.na(.), 0L, .))
    )
}


# -----------------------------------------------------------------------------
# 4. SUMMARY STATISTICS
# -----------------------------------------------------------------------------
n_cells    <- nrow(calendar)
n_with_any <- sum(calendar$has_ndvi | calendar$has_lst)
n_flagged  <- sum(calendar$status == "Flagged: low quality")
pct_cover  <- round(100 * n_with_any / n_cells, 1)

annual <- calendar |>
  dplyr::group_by(year) |>
  dplyr::summarise(
    months_complete = sum(has_ndvi & has_lst),
    months_partial  = sum(xor(has_ndvi, has_lst)),
    months_empty    = sum(!has_ndvi & !has_lst),
    months_flagged  = sum(status == "Flagged: low quality"),
    .groups = "drop"
  )

# A year whose available months are all clustered in one season yields a biased
# annual mean even when the month count looks adequate, so the seasonal spread
# is reported alongside the count. Season labels follow HEMISPHERE
# (00_common.R, section 1e); when that is NA the column is omitted rather than
# fabricated, since a four-season scheme does not describe every climate.
SEASONS_ENABLED <- length(season_levels()) > 0

if (SEASONS_ENABLED) {
  annual <- annual |>
    dplyr::left_join(
      calendar |>
        dplyr::filter(has_ndvi | has_lst) |>
        dplyr::mutate(season = .season_of(month)) |>
        dplyr::group_by(year) |>
        dplyr::summarise(seasons_represented = dplyr::n_distinct(season),
                         .groups = "drop"),
      by = "year"
    ) |>
    dplyr::mutate(seasons_represented = ifelse(is.na(seasons_represented), 0L,
                                               seasons_represented))
}


# -----------------------------------------------------------------------------
# 5. WRITE TABLES
# -----------------------------------------------------------------------------
readr::write_csv(
  dplyr::select(calendar, year, month, status, has_ndvi, has_lst,
                dplyr::any_of(c("n_scenes", "sensors_used"))),
  file.path(DIR_OUT, "coverage_calendar_city.csv"))

readr::write_csv(annual, file.path(DIR_OUT, "coverage_summary_annual.csv"))

if (!is.null(zone_cover))
  readr::write_csv(zone_cover,
                   file.path(DIR_OUT, "coverage_calendar_zones.csv"))


# -----------------------------------------------------------------------------
# 6. TEXT REPORT
# -----------------------------------------------------------------------------
report <- file.path(DIR_OUT, "coverage_report.txt")
open_report(report)

cat(strrep("=", 78), "\n")
cat("  NDVI / LST COVERAGE CALENDAR -", AOI_LABEL, "\n")
cat(sprintf("  Period: %d-%d  (%d years x 12 months = %d cells)\n",
            YEAR_START, YEAR_END, YEAR_END - YEAR_START + 1, n_cells))
cat(sprintf("  Generated: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(strrep("=", 78), "\n\n")

cat("--- OVERALL ---\n")
cat(sprintf("  Cells with data (NDVI or LST) : %d / %d  (%.1f%%)\n",
            n_with_any, n_cells, pct_cover))
cat(sprintf("  Cells flagged low quality     : %d\n", n_flagged))
cat(sprintf("  Cells with no data at all     : %d\n",
            sum(calendar$status == "No data")))
cat("\n")

cat("--- BY YEAR ---\n")
cat("  months_complete     = both NDVI and LST present\n")
if (SEASONS_ENABLED) {
  cat("  seasons_represented = how many of the four seasons contribute; a year\n")
  cat("                        with data in only one or two seasons yields a\n")
  cat("                        seasonally biased annual mean regardless of how\n")
  cat("                        many months it has.\n")
}
cat("\n")
print(as.data.frame(annual), row.names = FALSE)
cat("\n")

cat("--- CALENDAR (rows = years, columns = months) ---\n")
cat("  Legend:  O = NDVI + LST | n = NDVI only | l = LST only\n")
cat("           ! = flagged low quality | . = no data\n\n")

symbol_for <- function(s) {
  dplyr::case_when(
    s == "NDVI + LST"           ~ "O",
    s == "NDVI only"            ~ "n",
    s == "LST only"             ~ "l",
    s == "Flagged: low quality" ~ "!",
    TRUE                        ~ "."
  )
}

wide <- calendar |>
  dplyr::mutate(sym = symbol_for(status)) |>
  dplyr::select(year, month, sym) |>
  tidyr::pivot_wider(names_from = month, values_from = sym, names_prefix = "M")

cat(sprintf("  year  %s\n", paste(sprintf("%2s", month.abb), collapse = " ")))
cat("  ", strrep("-", 4 + 12 * 3), "\n", sep = "")
for (i in seq_len(nrow(wide))) {
  vals <- as.character(wide[i, paste0("M", 1:12)])
  cat(sprintf("  %d  %s\n", wide$year[i],
              paste(sprintf("%2s", vals), collapse = " ")))
}
cat("\n")

if (!is.null(zone_cover)) {
  cat("--- ZONE-LEVEL COVERAGE (% of zones with an LST value) ---\n")
  cat("  A month can hold a scene and still leave zones empty, because cloud\n")
  cat("  or Landsat-7 SLC-off stripes may cover only part of the city. Values\n")
  cat("  well below 100 mark months whose zone-level statistics rest on a\n")
  cat("  subset of the city.\n\n")
  cw <- zone_cover |>
    dplyr::select(year, month, pct_lst) |>
    tidyr::pivot_wider(names_from = month, values_from = pct_lst,
                       names_prefix = "M")
  cat(sprintf("  year  %s\n", paste(sprintf("%5s", month.abb), collapse = " ")))
  cat("  ", strrep("-", 4 + 12 * 6), "\n", sep = "")
  for (i in seq_len(nrow(cw))) {
    vals <- as.numeric(cw[i, paste0("M", 1:12)])
    cat(sprintf("  %d  %s\n", cw$year[i],
                paste(ifelse(is.na(vals), "   NA", sprintf("%5.1f", vals)),
                      collapse = " ")))
  }
  cat("\n")
}

cat("--- FLAGGED MONTHS ---\n")
cat(excluded_months_report())
cat("\n")
cat(strrep("=", 78), "\n")

close_report()
message(sprintf("%s Report: %s", TAG, basename(report)))


# -----------------------------------------------------------------------------
# 7. FIGURE 1 - CITY COVERAGE CALENDAR
# -----------------------------------------------------------------------------
# Design: a categorical heat map. Years run down the vertical axis with the
# earliest at the top, matching how a reader scans a timeline. The number
# inside each cell is the count of scenes that were composited into that month,
# which distinguishes a month resting on one scene from a well-sampled one.

cal_plot <- calendar |>
  dplyr::mutate(
    status = factor(status,
                    levels = c("NDVI + LST", "NDVI only", "LST only",
                               "Flagged: low quality", "No data")),
    month_abb = factor(month.abb[month], levels = month.abb)
  )

status_cols <- c(
  "NDVI + LST"           = "#1a7a3e",
  "NDVI only"            = "#94c47d",
  "LST only"             = "#f0a860",
  "Flagged: low quality" = "#c0392b",
  "No data"              = "grey90"
)

has_scene_col <- "n_scenes" %in% names(cal_plot)

p1 <- ggplot(cal_plot, aes(x = month_abb, y = factor(year), fill = status)) +
  geom_tile(colour = "white", linewidth = 0.8)

if (has_scene_col) {
  p1 <- p1 + geom_text(
    aes(label = ifelse(!is.na(n_scenes) & n_scenes > 0,
                       as.character(n_scenes), ""),
        # White text on the dark green and dark red fills, dark text elsewhere:
        # a single text colour would be unreadable on one end of the scale.
        colour = status %in% c("NDVI + LST", "Flagged: low quality")),
    size = 3, fontface = "bold", show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20"))
}

p1 <- p1 +
  scale_fill_manual(values = status_cols, name = NULL, drop = FALSE) +
  scale_y_discrete(limits = rev) +
  coord_equal() +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                             override.aes = list(colour = NA))) +
  labs(
    title    = sprintf("Monthly availability of NDVI and LST, %s", AOI_LABEL),
    subtitle = sprintf(
      "%d of %d month-cells carry data (%.0f%%). Numbers give the count of Landsat scenes composited into each month.",
      n_with_any, n_cells, pct_cover),
    x = NULL, y = "Year",
    caption = paste0(
      "Grey cells are months with no cloud-free acquisition. ",
      if (nzchar(excluded_months_caption()))
        paste0("Red marks ", excluded_months_caption(),
               "from all downstream statistics ")
      else "",
      "(Landsat-7 SLC-off artefact).\n",
      "Data: Landsat Collection 2 Level-2 (USGS).")
  ) +
  theme_pub(grid = "none") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = FIG_BASE_SIZE - 2))

save_fig(p1, "fig01_coverage_calendar_city", DIR_OUT,
         width = 11, height = 8.5, tag = TAG)


# -----------------------------------------------------------------------------
# 8. FIGURE 2 - ZONE-LEVEL SPATIAL COVERAGE
# -----------------------------------------------------------------------------
# Design: a continuous heat map of the percentage of zones with an LST value.
# The sequential single-hue ramp is deliberate - this is a magnitude, not a
# diverging quantity, so a diverging palette would imply a meaningless midpoint.

if (!is.null(zone_cover)) {

  zc_plot <- zone_cover |>
    dplyr::mutate(month_abb = factor(month.abb[month], levels = month.abb))

  p2 <- ggplot(zc_plot, aes(x = month_abb, y = factor(year), fill = pct_lst)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = ifelse(pct_lst > 0, sprintf("%.0f", pct_lst), ""),
                  colour = pct_lst > 55),
              size = 2.7, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey25")) +
    scale_fill_gradientn(
      colours = c("grey92", "#deebf7", "#9ecae1", "#4292c6", "#08519c"),
      limits  = c(0, 100),
      breaks  = c(0, 25, 50, 75, 100),
      labels  = function(x) paste0(x, "%"),
      name    = "Zones with\nLST data"
    ) +
    scale_y_discrete(limits = rev) +
    coord_equal() +
    labs(
      title    = sprintf("Spatial completeness of the LST record, %s", AOI_LABEL),
      subtitle = sprintf(
        "Percentage of the %d zones holding at least one valid LST pixel in each month.",
        max(zone_cover$n_zones_total, na.rm = TRUE)),
      x = NULL, y = "Year",
      caption = paste0(
        "Values below 100% mark months in which a scene exists but cloud or ",
        "Landsat-7 SLC-off stripes leave part of the city unobserved;\n",
        "zone-level statistics for those months rest on a subset of the study area.")
    ) +
    theme_pub(grid = "none") +
    theme(legend.position = "right",
          axis.text.x = element_text(size = FIG_BASE_SIZE - 2))

  save_fig(p2, "fig02_coverage_calendar_zones", DIR_OUT,
           width = 11, height = 8.5, tag = TAG)
}


# -----------------------------------------------------------------------------
# 9. FIGURE 3 - ANNUAL DATA AVAILABILITY
# -----------------------------------------------------------------------------
# A stacked bar per year, showing how the twelve months break down. This is the
# figure to point at when justifying why some years are excluded, or why an
# annual mean carries a caveat.

annual_long <- annual |>
  dplyr::select(year, months_complete, months_partial, months_empty) |>
  tidyr::pivot_longer(-year, names_to = "category", values_to = "n_months") |>
  dplyr::mutate(category = factor(
    dplyr::recode(category,
                  months_complete = "NDVI + LST",
                  months_partial  = "One variable only",
                  months_empty    = "No data"),
    levels = c("NDVI + LST", "One variable only", "No data")))

p3 <- ggplot(annual_long, aes(x = factor(year), y = n_months, fill = category)) +
  geom_col(width = 0.78) +
  # A reference line at six months: below this the annual mean rests on half a
  # year or less and is unlikely to be seasonally balanced.
  geom_hline(yintercept = 6, linetype = "dashed",
             colour = "grey30", linewidth = 0.4) +
  annotate("text", x = 0.7, y = 6.45, label = "6 months",
           hjust = 0, size = 2.9, colour = "grey30") +
  scale_fill_manual(values = c("NDVI + LST"        = "#1a7a3e",
                               "One variable only" = "#f0a860",
                               "No data"           = "grey88"),
                    name = NULL) +
  scale_y_continuous(breaks = seq(0, 12, 2), expand = expansion(mult = c(0, 0.04))) +
  labs(
    title    = sprintf("Annual data availability, %s", AOI_LABEL),
    subtitle = "Composition of each year's twelve months. Annual means from years below the dashed line should be read with caution.",
    x = NULL, y = "Number of months",
    caption = if (SEASONS_ENABLED)
      "Seasonal balance matters as much as the count: see coverage_summary_annual.csv for the number of seasons represented in each year."
    else
      "See coverage_summary_annual.csv for the month-by-month breakdown."
  ) +
  theme_pub(grid = "y") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

save_fig(p3, "fig03_annual_availability", DIR_OUT,
         width = 10, height = 5.5, tag = TAG)


message(sprintf("%s Done: %s", TAG, Sys.time()))
message(sprintf("%s Outputs in %s\n", TAG, DIR_OUT))
