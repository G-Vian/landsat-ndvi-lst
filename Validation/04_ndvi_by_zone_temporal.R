# =============================================================================
# VALIDATION - 04 - SPATIO-TEMPORAL VARIATION OF NDVI ACROSS ZONES
# =============================================================================
# QUESTION ANSWERED
#   "How does vegetation cover vary across the zones of the city, and within
#    any one zone, is it stable through time or does it change?"
#
# This script mirrors 03_lst_by_zone_temporal.R, with one deliberate
# methodological difference explained below.
#
# INTERPRETING NDVI VALUES
#   NDVI is bounded to [-1, 1] and is dimensionless. Conventional bands:
#     below 0     water, and occasionally wet asphalt or deep shadow
#     0 to 0.2    bare soil, sand, dense built-up fabric with no canopy
#     0.2 to 0.4  sparse vegetation, dry grass, mixed urban
#     0.4 to 0.6  moderate vegetation, tree-lined streets, parks
#     above 0.6   dense vegetation, closed forest canopy
#   In a coastal city such as Santos most zones fall between 0.2 and 0.5.
#
# WHY NO COEFFICIENT OF VARIATION HERE
#   The LST script reports a coefficient of variation, sd divided by mean.
#   That statistic is only meaningful when the mean is comfortably away from
#   zero. NDVI means can legitimately approach zero in heavily built zones, at
#   which point the ratio explodes and produces enormous values that reflect a
#   near-zero denominator rather than any real instability. This script
#   therefore reports the plain standard deviation and the observed range,
#   which stay interpretable across the whole NDVI scale.
#
# HOW TO RUN
#   singularity exec --bind /home/g.vian:/home/g.vian \
#     /home/public/R_inla/r_inla.sif \
#     Rscript 04_ndvi_by_zone_temporal.R
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

TAG     <- "[04_ndvi_zones]"
DIR_OUT <- file.path(VALIDATION_ROOT, "04_ndvi_by_zone")
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

MIN_MONTHS_FOR_TREND <- 24

message(sprintf("\n%s Start: %s", TAG, Sys.time()))
message(sprintf("%s Output: %s", TAG, DIR_OUT))


# -----------------------------------------------------------------------------
# 1. LOAD AND CLEAN
# -----------------------------------------------------------------------------
pipeline_dir <- resolve_pipeline_dir()
tables_dir   <- resolve_tables_dir(pipeline_dir)

df <- load_zone_monthly(tables_dir) |>
  dplyr::select(dplyr::any_of(c("zone_id", "zone_name", "year", "month",
                                "ndvi_mean", "ndvi_sd", "ndvi_min",
                                "ndvi_max", "ndvi_n", "n_scenes"))) |>
  add_time_columns()

df_clean <- df |> dplyr::filter(!is.na(ndvi_mean), !low_quality)

n_zones  <- dplyr::n_distinct(df_clean$zone_name)
n_months <- dplyr::n_distinct(df_clean$date)
message(sprintf("%s %d zones, %d months with data.", TAG, n_zones, n_months))


# -----------------------------------------------------------------------------
# 2. PER-ZONE STATISTICS
# -----------------------------------------------------------------------------
# ORDER OF OPERATIONS MATTERS INSIDE summarise().
# Every statistic that needs the VECTOR of monthly values - median, sd, min,
# max, the trend - has to be computed before any line that reassigns the name
# ndvi_mean. From the moment a summarise() line writes to ndvi_mean, that name
# refers to a single number rather than to the column, and every later
# reference silently returns that scalar. Here the reassignment is avoided
# altogether by giving the zone-level mean a distinct name.

zone_stats <- df_clean |>
  dplyr::group_by(zone_id, zone_name) |>
  dplyr::summarise(
    n_months_obs  = dplyr::n(),
    ndvi_mean_all = mean(ndvi_mean, na.rm = TRUE),
    ndvi_median   = stats::median(ndvi_mean, na.rm = TRUE),
    ndvi_sd_time  = stats::sd(ndvi_mean, na.rm = TRUE),
    ndvi_min_all  = min(ndvi_mean, na.rm = TRUE),
    ndvi_max_all  = max(ndvi_mean, na.rm = TRUE),
    ndvi_range    = max(ndvi_mean, na.rm = TRUE) - min(ndvi_mean, na.rm = TRUE),
    trend_per_year = slope_per_year(date, ndvi_mean,
                                    min_n = MIN_MONTHS_FOR_TREND),
    .groups = "drop"
  ) |>
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~round(., 4))) |>
  dplyr::arrange(dplyr::desc(ndvi_mean_all))

readr::write_csv(zone_stats, file.path(DIR_OUT, "statistics_by_zone.csv"))

climatology <- df_clean |>
  dplyr::group_by(zone_name, month) |>
  dplyr::summarise(ndvi_month_mean = round(mean(ndvi_mean, na.rm = TRUE), 4),
                   n_obs = dplyr::n(), .groups = "drop")
readr::write_csv(climatology, file.path(DIR_OUT, "climatology_by_zone.csv"))

annual <- df_clean |>
  dplyr::group_by(zone_name, year) |>
  dplyr::summarise(ndvi_year_mean = mean(ndvi_mean, na.rm = TRUE),
                   n_obs = dplyr::n(), .groups = "drop")
readr::write_csv(annual, file.path(DIR_OUT, "annual_means_by_zone.csv"))


# -----------------------------------------------------------------------------
# 3. TEXT REPORT
# -----------------------------------------------------------------------------
greenest  <- utils::head(dplyr::arrange(zone_stats, dplyr::desc(ndvi_mean_all)), 5)
barest    <- utils::head(dplyr::arrange(zone_stats, ndvi_mean_all), 5)
most_var  <- utils::head(dplyr::arrange(zone_stats, dplyr::desc(ndvi_sd_time)), 5)
least_var <- utils::head(dplyr::arrange(zone_stats, ndvi_sd_time), 5)

trend_ok  <- zone_stats |> dplyr::filter(!is.na(trend_per_year))
losing    <- utils::head(dplyr::arrange(trend_ok, trend_per_year), 5)
gaining   <- utils::head(dplyr::arrange(trend_ok, dplyr::desc(trend_per_year)), 5)

report <- file.path(DIR_OUT, "ndvi_by_zone_report.txt")
open_report(report)

cat(strrep("=", 78), "\n")
cat("  SPATIO-TEMPORAL VARIATION OF NDVI -", AOI_LABEL, "\n")
cat(sprintf("  Zones: %d | Months: %d | Generated: %s\n",
            n_zones, n_months, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("  Excluded: %d-%02d\n", LOW_QUALITY_MONTHS$year[1],
            LOW_QUALITY_MONTHS$month[1]))
cat(strrep("=", 78), "\n\n")

cat("--- CITY-WIDE SUMMARY ---\n")
cat(sprintf("  Mean NDVI across zones : %.3f\n",
            mean(zone_stats$ndvi_mean_all)))
cat(sprintf("  Greenest zone : %s (NDVI %.3f)\n",
            greenest$zone_name[1], greenest$ndvi_mean_all[1]))
cat(sprintf("  Barest zone   : %s (NDVI %.3f)\n",
            barest$zone_name[1], barest$ndvi_mean_all[1]))
cat(sprintf("  Spread between the two : %.3f\n",
            greenest$ndvi_mean_all[1] - barest$ndvi_mean_all[1]))
cat("\n")

cat("--- REFERENCE VALUES FOR NDVI ---\n")
cat("  below 0    : water, wet asphalt, deep shadow\n")
cat("  0.0 - 0.2  : bare soil, sand, dense built-up with no canopy\n")
cat("  0.2 - 0.4  : sparse vegetation, dry grass, mixed urban\n")
cat("  0.4 - 0.6  : moderate vegetation, tree-lined streets, parks\n")
cat("  above 0.6  : dense vegetation, closed canopy\n\n")

cat("--- FIVE MOST VEGETATED ZONES ---\n")
print(as.data.frame(greenest[, c("zone_name", "ndvi_mean_all", "ndvi_sd_time",
                                 "ndvi_range", "n_months_obs")]),
      row.names = FALSE)
cat("\n--- FIVE LEAST VEGETATED ZONES ---\n")
print(as.data.frame(barest[, c("zone_name", "ndvi_mean_all", "ndvi_sd_time",
                               "ndvi_range", "n_months_obs")]),
      row.names = FALSE)
cat("\n")

cat("--- MOST VARIABLE THROUGH TIME ---\n")
print(as.data.frame(most_var[, c("zone_name", "ndvi_mean_all", "ndvi_sd_time",
                                 "ndvi_min_all", "ndvi_max_all")]),
      row.names = FALSE)
cat("\n--- MOST STABLE THROUGH TIME ---\n")
print(as.data.frame(least_var[, c("zone_name", "ndvi_mean_all", "ndvi_sd_time",
                                  "ndvi_min_all", "ndvi_max_all")]),
      row.names = FALSE)
cat("\n")

cat("--- LONG-TERM TREND ---\n")
cat("  Slope of a regression on decimal year, in NDVI units per year.\n")
cat("  Positive means greening; negative means loss of vegetation, typically\n")
cat("  through construction or canopy removal.\n\n")
if (nrow(trend_ok) > 0) {
  cat(sprintf("  Median trend across zones : %+.4f per year\n",
              stats::median(trend_ok$trend_per_year, na.rm = TRUE)))
  cat(sprintf("  Strongest greening : %+.4f per year (%s)\n",
              max(trend_ok$trend_per_year, na.rm = TRUE),
              trend_ok$zone_name[which.max(trend_ok$trend_per_year)]))
  cat(sprintf("  Strongest loss     : %+.4f per year (%s)\n\n",
              min(trend_ok$trend_per_year, na.rm = TRUE),
              trend_ok$zone_name[which.min(trend_ok$trend_per_year)]))
  cat("  Five zones losing vegetation fastest:\n")
  print(as.data.frame(losing[, c("zone_name", "trend_per_year",
                                 "ndvi_mean_all", "n_months_obs")]),
        row.names = FALSE)
  cat("\n  Five zones gaining vegetation fastest:\n")
  print(as.data.frame(gaining[, c("zone_name", "trend_per_year",
                                  "ndvi_mean_all", "n_months_obs")]),
        row.names = FALSE)
} else {
  cat("  No zone has enough observations for a stable trend estimate.\n")
}
cat("\n")

cat("--- CAVEATS ---\n")
cat(sprintf("  Trends use only zones with at least %d monthly observations.\n",
            MIN_MONTHS_FOR_TREND))
cat("  NDVI has a marked seasonal cycle, so a trend estimated from a record\n")
cat("  whose available months are unevenly distributed across the year can\n")
cat("  reflect that imbalance rather than any real change in vegetation.\n")
cat("  Cross-check against the coverage calendar before interpreting a slope.\n\n")
cat(strrep("=", 78), "\n")

close_report()
message(sprintf("%s Report: %s", TAG, basename(report)))


# -----------------------------------------------------------------------------
# 4. FIGURES
# -----------------------------------------------------------------------------

city_mean <- df_clean |>
  dplyr::group_by(date) |>
  dplyr::summarise(ndvi_city = mean(ndvi_mean, na.rm = TRUE), .groups = "drop")

# --- Figure 1: time series, one line per zone --------------------------------
# As in the LST script, lines are coloured by the zone's long-run mean so that
# the palette encodes a real variable rather than arbitrary zone identity.

spaghetti <- df_clean |>
  dplyr::left_join(dplyr::select(zone_stats, zone_name, ndvi_mean_all),
                   by = "zone_name")

p1 <- ggplot() +
  geom_line(data = spaghetti,
            aes(x = date, y = ndvi_mean, group = zone_name,
                colour = ndvi_mean_all),
            linewidth = 0.32, alpha = 0.55) +
  geom_line(data = city_mean, aes(x = date, y = ndvi_city),
            colour = "grey12", linewidth = 1.1) +
  geom_point(data = city_mean, aes(x = date, y = ndvi_city),
             colour = "grey12", size = 1.1) +
  scale_colour_gradientn(colours = PAL_NDVI, name = "Zone long-run\nmean NDVI") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y",
               expand = expansion(mult = 0.01)) +
  labs(
    title    = sprintf("Monthly NDVI by zone, %s", AOI_LABEL),
    subtitle = sprintf("%d zones, coloured by their long-run mean. The heavy dark line is the city-wide average.",
                       n_zones),
    x = NULL, y = "NDVI (dimensionless)",
    caption = sprintf(
      "Green lines stay in the upper part of the bundle throughout, indicating that the spatial pattern of vegetation is stable over the record.\nJuly 2021 excluded. Data: Landsat Collection 2 Level-2 (USGS), %d-%d.",
      YEAR_START, YEAR_END)
  ) +
  theme_pub(grid = "y") +
  theme(legend.position = "right")

save_fig(p1, "fig01_timeseries_by_zone", DIR_OUT, 12, 6, TAG)


# --- Figure 2: distribution per zone -----------------------------------------

box_df <- df_clean |>
  dplyr::group_by(zone_name) |>
  dplyr::mutate(zone_median = stats::median(ndvi_mean, na.rm = TRUE)) |>
  dplyr::ungroup() |>
  dplyr::mutate(zone_name = stats::reorder(zone_name, zone_median))

p2 <- ggplot(box_df, aes(x = ndvi_mean, y = zone_name, fill = zone_median)) +
  geom_boxplot(outlier.size = 0.55, outlier.alpha = 0.45, linewidth = 0.25,
               width = 0.7) +
  scale_fill_gradientn(colours = PAL_NDVI, name = "Median\nNDVI") +
  labs(
    title    = sprintf("Distribution of monthly NDVI by zone, %s", AOI_LABEL),
    subtitle = "Zones ordered by median, from the most built-up at the bottom to the most vegetated at the top.",
    x = "NDVI (dimensionless)", y = NULL,
    caption = "Each box summarises every month in which the zone was observed."
  ) +
  theme_pub(base_size = 9, grid = "x") +
  theme(legend.position = "right")

save_fig(p2, "fig02_distribution_by_zone", DIR_OUT, 9, 12.5, TAG)


# --- Figure 3: monthly climatology heat map ----------------------------------

p3 <- climatology |>
  dplyr::mutate(month_abb = factor(month.abb[month], levels = month.abb),
                zone_name = stats::reorder(zone_name, ndvi_month_mean,
                                           FUN = mean)) |>
  ggplot(aes(x = month_abb, y = zone_name, fill = ndvi_month_mean)) +
  geom_tile(colour = "white", linewidth = 0.12) +
  scale_fill_gradientn(colours = PAL_NDVI, name = "Mean\nNDVI") +
  scale_x_discrete(expand = c(0, 0)) +
  labs(
    title    = sprintf("Seasonal cycle of NDVI by zone, %s", AOI_LABEL),
    subtitle = "Mean value of each zone in each calendar month, pooling all years.",
    x = NULL, y = NULL,
    caption = "Strong horizontal banding indicates that differences between zones are much larger than the seasonal variation within any one zone."
  ) +
  theme_pub(base_size = 9, grid = "none") +
  theme(legend.position = "right")

save_fig(p3, "fig03_climatology_heatmap", DIR_OUT, 9, 12.5, TAG)


# --- Figure 4: interannual heat map ------------------------------------------

p4 <- annual |>
  dplyr::mutate(zone_name = stats::reorder(zone_name, ndvi_year_mean,
                                           FUN = mean)) |>
  ggplot(aes(x = factor(year), y = zone_name, fill = ndvi_year_mean)) +
  geom_tile(colour = "white", linewidth = 0.12) +
  scale_fill_gradientn(colours = PAL_NDVI, name = "Mean\nNDVI",
                       na.value = "grey92") +
  scale_x_discrete(expand = c(0, 0)) +
  labs(
    title    = sprintf("Annual mean NDVI by zone, %s", AOI_LABEL),
    subtitle = "Blank cells are zone-years with no valid observation.",
    x = NULL, y = NULL,
    caption = "A zone whose row changes colour progressively from left to right is gaining or losing vegetation; see the trend figure for the fitted rate."
  ) +
  theme_pub(base_size = 9, grid = "none") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))

save_fig(p4, "fig04_annual_heatmap", DIR_OUT, 9.5, 12.5, TAG)


# --- Figure 5: temporal variability by zone ----------------------------------

p5 <- zone_stats |>
  dplyr::arrange(ndvi_sd_time) |>
  dplyr::mutate(zone_name = factor(zone_name, levels = zone_name)) |>
  ggplot(aes(x = ndvi_sd_time, y = zone_name, fill = ndvi_sd_time)) +
  geom_col(width = 0.75) +
  scale_fill_gradient(low = "#66bd63", high = "#a50026",
                      name = "sd of NDVI") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.03))) +
  labs(
    title    = sprintf("Temporal variability of NDVI by zone, %s", AOI_LABEL),
    subtitle = "Standard deviation of the monthly series within each zone. The plain standard deviation is used rather than a coefficient of variation, which is unstable when the mean approaches zero.",
    x = "Standard deviation across months (NDVI units)", y = NULL,
    caption = "Zones with a mixed surface of vegetation and built fabric vary most, because the vegetated fraction responds to season while the built fraction does not."
  ) +
  theme_pub(base_size = 9, grid = "x") +
  theme(legend.position = "right")

save_fig(p5, "fig05_variability_by_zone", DIR_OUT, 8.5, 12.5, TAG)


# --- Figure 6: long-term trend by zone ---------------------------------------

if (nrow(trend_ok) > 0) {
  p6 <- trend_ok |>
    dplyr::arrange(trend_per_year) |>
    dplyr::mutate(zone_name = factor(zone_name, levels = zone_name)) |>
    ggplot(aes(x = trend_per_year, y = zone_name, fill = trend_per_year)) +
    geom_col(width = 0.75) +
    geom_vline(xintercept = 0, colour = "grey25", linewidth = 0.4) +
    scale_fill_gradient2(low = "#a50026", mid = "grey95", high = "#1a9850",
                         midpoint = 0, name = "NDVI / year") +
    labs(
      title    = sprintf("Linear trend in NDVI by zone, %s (%d-%d)",
                         AOI_LABEL, YEAR_START, YEAR_END),
      subtitle = "Positive values indicate greening; negative values indicate loss of vegetation cover.",
      x = "Trend (NDVI units per year)", y = NULL,
      caption = sprintf(
        "Zones with fewer than %d monthly observations are omitted. These are descriptive slopes with no significance testing; NDVI is strongly seasonal, so uneven month coverage within a year can masquerade as a trend.",
        MIN_MONTHS_FOR_TREND)
    ) +
    theme_pub(base_size = 9, grid = "x") +
    theme(legend.position = "right")

  save_fig(p6, "fig06_trend_by_zone", DIR_OUT, 8.5, 12.5, TAG)
}


message(sprintf("%s Done: %s", TAG, Sys.time()))
message(sprintf("%s Outputs in %s\n", TAG, DIR_OUT))
