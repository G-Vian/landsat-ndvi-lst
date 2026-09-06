# =============================================================================
# VALIDATION - 03 - SPATIO-TEMPORAL VARIATION OF LST ACROSS ZONES
# =============================================================================
# QUESTION ANSWERED
#   "How does surface temperature vary across the zones of the city, and within
#    any one zone, is it roughly constant through time or does it swing?"
#
# WHY THIS MATTERS
#   A spatial model treats each zone as an exposure unit, so two properties of
#   the covariate have to be established before it is used. First, that zones
#   genuinely differ from one another - if every zone had the same surface
#   temperature there would be no spatial signal to exploit. Second, that the
#   within-zone temporal variation is understood, because a zone whose value
#   swings wildly from month to month contributes a noisy exposure even if its
#   long-run mean is well estimated.
#
# WHAT THE SCRIPT COMPUTES, PER ZONE
#   mean, median, standard deviation, minimum, maximum
#   coefficient of variation (sd / mean x 100), a scale-free measure of how
#     large the temporal swing is relative to the zone's own level
#   linear trend in degrees per year, from a regression on decimal year
#   monthly climatology, that is the average value in each calendar month
#
# HOW TO RUN
#   singularity exec --bind /home/g.vian:/home/g.vian \
#     /home/public/R_inla/r_inla.sif \
#     Rscript 03_lst_by_zone_temporal.R
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

TAG     <- "[03_lst_zones]"
DIR_OUT <- file.path(VALIDATION_ROOT, "03_lst_by_zone")
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

# A trend fitted to a short record is dominated by which months happen to be
# present rather than by any real change, so zones below this threshold get a
# missing slope instead of a misleading number.
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
                                "lst_mean", "lst_sd", "lst_min", "lst_max",
                                "lst_n", "n_scenes"))) |>
  add_time_columns()

df_clean <- df |> dplyr::filter(!is.na(lst_mean), !low_quality)

n_zones  <- dplyr::n_distinct(df_clean$zone_name)
n_months <- dplyr::n_distinct(df_clean$date)
message(sprintf("%s %d zones, %d months with data.", TAG, n_zones, n_months))


# -----------------------------------------------------------------------------
# 2. PER-ZONE STATISTICS
# -----------------------------------------------------------------------------
# Note on the coefficient of variation: it is only meaningful because surface
# temperature in Celsius is comfortably away from zero here. For a variable
# that can approach zero the ratio explodes and the plain standard deviation
# should be used instead - which is exactly the situation in the NDVI script.

zone_stats <- df_clean |>
  dplyr::group_by(zone_id, zone_name) |>
  dplyr::summarise(
    n_months_obs = dplyr::n(),
    lst_mean_all = mean(lst_mean, na.rm = TRUE),
    lst_median   = stats::median(lst_mean, na.rm = TRUE),
    lst_sd_time  = stats::sd(lst_mean, na.rm = TRUE),
    lst_min_all  = min(lst_mean, na.rm = TRUE),
    lst_max_all  = max(lst_mean, na.rm = TRUE),
    lst_range    = max(lst_mean, na.rm = TRUE) - min(lst_mean, na.rm = TRUE),
    lst_cv_pct   = 100 * stats::sd(lst_mean, na.rm = TRUE) /
                         mean(lst_mean, na.rm = TRUE),
    trend_c_per_year = slope_per_year(date, lst_mean,
                                      min_n = MIN_MONTHS_FOR_TREND),
    .groups = "drop"
  ) |>
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~round(., 3))) |>
  dplyr::arrange(dplyr::desc(lst_mean_all))

readr::write_csv(zone_stats, file.path(DIR_OUT, "statistics_by_zone.csv"))

# Monthly climatology: the mean of each zone in each calendar month, pooling
# all years. This separates the seasonal component from the interannual one.
climatology <- df_clean |>
  dplyr::group_by(zone_name, month) |>
  dplyr::summarise(lst_month_mean = round(mean(lst_mean, na.rm = TRUE), 2),
                   n_obs = dplyr::n(), .groups = "drop")
readr::write_csv(climatology, file.path(DIR_OUT, "climatology_by_zone.csv"))

# Annual means per zone, used for the interannual heat map.
annual <- df_clean |>
  dplyr::group_by(zone_name, year) |>
  dplyr::summarise(lst_year_mean = mean(lst_mean, na.rm = TRUE),
                   n_obs = dplyr::n(), .groups = "drop")
readr::write_csv(annual, file.path(DIR_OUT, "annual_means_by_zone.csv"))


# -----------------------------------------------------------------------------
# 3. TEXT REPORT
# -----------------------------------------------------------------------------
hottest   <- utils::head(dplyr::arrange(zone_stats, dplyr::desc(lst_mean_all)), 5)
coolest   <- utils::head(dplyr::arrange(zone_stats, lst_mean_all), 5)
most_var  <- utils::head(dplyr::arrange(zone_stats, dplyr::desc(lst_sd_time)), 5)
least_var <- utils::head(dplyr::arrange(zone_stats, lst_sd_time), 5)
warming   <- zone_stats |> dplyr::filter(!is.na(trend_c_per_year)) |>
  dplyr::arrange(dplyr::desc(trend_c_per_year)) |> utils::head(5)

report <- file.path(DIR_OUT, "lst_by_zone_report.txt")
open_report(report)

cat(strrep("=", 78), "\n")
cat("  SPATIO-TEMPORAL VARIATION OF SURFACE TEMPERATURE -", AOI_LABEL, "\n")
cat(sprintf("  Zones: %d | Months: %d | Generated: %s\n",
            n_zones, n_months, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("  Excluded: %d-%02d (%s)\n", LOW_QUALITY_MONTHS$year[1],
            LOW_QUALITY_MONTHS$month[1], LOW_QUALITY_MONTHS$reason[1]))
cat(strrep("=", 78), "\n\n")

cat("--- CITY-WIDE SUMMARY ---\n")
cat(sprintf("  Mean surface temperature across zones : %.2f C\n",
            mean(zone_stats$lst_mean_all)))
cat(sprintf("  Warmest zone : %s (%.2f C)\n",
            hottest$zone_name[1], hottest$lst_mean_all[1]))
cat(sprintf("  Coolest zone : %s (%.2f C)\n",
            coolest$zone_name[1], coolest$lst_mean_all[1]))
cat(sprintf("  Spread between the two : %.2f C\n",
            hottest$lst_mean_all[1] - coolest$lst_mean_all[1]))
cat("\n")
cat("  This spread is the spatial contrast a zone-level model exploits. It\n")
cat("  should be compared against the measurement error reported by the\n")
cat("  inter-sensor analysis: the contrast has to be comfortably larger than\n")
cat("  the error for zone-level effects to be estimable.\n\n")

cat("--- FIVE WARMEST ZONES (long-run mean) ---\n")
print(as.data.frame(hottest[, c("zone_name", "lst_mean_all", "lst_sd_time",
                                "lst_cv_pct", "n_months_obs")]),
      row.names = FALSE)
cat("\n--- FIVE COOLEST ZONES ---\n")
print(as.data.frame(coolest[, c("zone_name", "lst_mean_all", "lst_sd_time",
                                "lst_cv_pct", "n_months_obs")]),
      row.names = FALSE)
cat("\n")

cat("--- MOST VARIABLE THROUGH TIME (largest temporal sd) ---\n")
print(as.data.frame(most_var[, c("zone_name", "lst_mean_all", "lst_sd_time",
                                 "lst_min_all", "lst_max_all")]),
      row.names = FALSE)
cat("\n--- MOST STABLE THROUGH TIME (smallest temporal sd) ---\n")
print(as.data.frame(least_var[, c("zone_name", "lst_mean_all", "lst_sd_time",
                                  "lst_min_all", "lst_max_all")]),
      row.names = FALSE)
cat("\n")

cat("--- LONG-TERM TREND ---\n")
cat(sprintf("  Fitted as a linear regression on decimal year, reported in\n"))
cat(sprintf("  degrees Celsius per year. Zones with fewer than %d monthly\n",
            MIN_MONTHS_FOR_TREND))
cat("  observations are left blank rather than given an unstable slope.\n\n")
if (nrow(warming) > 0) {
  cat(sprintf("  Median trend across zones : %+.3f C/year\n",
              stats::median(zone_stats$trend_c_per_year, na.rm = TRUE)))
  cat(sprintf("  Strongest warming : %+.3f C/year (%s)\n",
              max(zone_stats$trend_c_per_year, na.rm = TRUE),
              zone_stats$zone_name[which.max(zone_stats$trend_c_per_year)]))
  cat(sprintf("  Strongest cooling : %+.3f C/year (%s)\n\n",
              min(zone_stats$trend_c_per_year, na.rm = TRUE),
              zone_stats$zone_name[which.min(zone_stats$trend_c_per_year)]))
  cat("  Five zones with the fastest warming:\n")
  print(as.data.frame(warming[, c("zone_name", "trend_c_per_year",
                                  "lst_mean_all", "n_months_obs")]),
        row.names = FALSE)
} else {
  cat("  No zone has enough observations for a stable trend estimate.\n")
}
cat("\n")

cat("--- HOW TO READ THE COEFFICIENT OF VARIATION ---\n")
cat("  CV = sd / mean x 100. It expresses the temporal swing relative to the\n")
cat("  zone's own level, so zones with different baselines can be compared.\n")
cat("  A high CV suggests a mixed surface whose apparent temperature responds\n")
cat("  strongly to weather; a low CV suggests a uniform surface, either dense\n")
cat("  vegetation or homogeneous built-up fabric.\n")
cat(sprintf("  Median CV across zones : %.2f%%\n\n",
            stats::median(zone_stats$lst_cv_pct, na.rm = TRUE)))

cat("--- CAVEAT ---\n")
cat("  Every statistic above is computed only from months in which the zone\n")
cat("  was observed. Zones differ in how many months that is, so a zone with\n")
cat("  a short record may appear unusually stable simply because it has been\n")
cat("  sampled in fewer distinct conditions. Always read n_months_obs\n")
cat("  alongside the statistics.\n\n")
cat(strrep("=", 78), "\n")

close_report()
message(sprintf("%s Report: %s", TAG, basename(report)))


# -----------------------------------------------------------------------------
# 4. FIGURES
# -----------------------------------------------------------------------------

city_mean <- df_clean |>
  dplyr::group_by(date) |>
  dplyr::summarise(lst_city = mean(lst_mean, na.rm = TRUE), .groups = "drop")

# --- Figure 1: time series, one line per zone --------------------------------
# The original version of this figure coloured the lines by zone identity,
# which with 60 unlabelled categories carries no information at all. Colouring
# instead by each zone's long-run mean turns the palette into an actual
# variable: the reader can see that warm zones stay near the top of the bundle
# throughout the record, which is the persistence result made visual.

spaghetti <- df_clean |>
  dplyr::left_join(dplyr::select(zone_stats, zone_name, lst_mean_all),
                   by = "zone_name")

p1 <- ggplot() +
  geom_line(data = spaghetti,
            aes(x = date, y = lst_mean, group = zone_name,
                colour = lst_mean_all),
            linewidth = 0.32, alpha = 0.55) +
  geom_line(data = city_mean, aes(x = date, y = lst_city),
            colour = "grey12", linewidth = 1.1) +
  geom_point(data = city_mean, aes(x = date, y = lst_city),
             colour = "grey12", size = 1.1) +
  scale_colour_gradientn(colours = PAL_LST, name = "Zone long-run\nmean (\u00b0C)") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y",
               expand = expansion(mult = 0.01)) +
  labs(
    title    = sprintf("Monthly surface temperature by zone, %s", AOI_LABEL),
    subtitle = sprintf("%d zones, coloured by their long-run mean. The heavy dark line is the city-wide average.",
                       n_zones),
    x = NULL, y = "Land surface temperature (\u00b0C)",
    caption = sprintf(
      "Warm-coloured lines remain in the upper part of the bundle throughout, indicating that the spatial ordering of zones persists across seasons.\nJuly 2021 excluded. Data: Landsat Collection 2 Level-2 (USGS), %d-%d.",
      YEAR_START, YEAR_END)
  ) +
  theme_pub(grid = "y") +
  theme(legend.position = "right")

save_fig(p1, "fig01_timeseries_by_zone", DIR_OUT, 12, 6, TAG)


# --- Figure 2: distribution per zone -----------------------------------------
# Ordered by median so the vertical axis itself carries the ranking. The fill
# repeats the median, which makes the warm/cool gradient legible at a glance
# even where the boxes are thin.

box_df <- df_clean |>
  dplyr::group_by(zone_name) |>
  dplyr::mutate(zone_median = stats::median(lst_mean, na.rm = TRUE)) |>
  dplyr::ungroup() |>
  dplyr::mutate(zone_name = stats::reorder(zone_name, zone_median))

p2 <- ggplot(box_df, aes(x = lst_mean, y = zone_name, fill = zone_median)) +
  geom_boxplot(outlier.size = 0.55, outlier.alpha = 0.45, linewidth = 0.25,
               width = 0.7) +
  scale_fill_gradientn(colours = PAL_LST, name = "Median\n(\u00b0C)") +
  labs(
    title    = sprintf("Distribution of monthly surface temperature by zone, %s", AOI_LABEL),
    subtitle = "Zones ordered by median. Box width shows the interquartile range of the monthly series; whiskers extend to 1.5 times that range.",
    x = "Land surface temperature (\u00b0C)", y = NULL,
    caption = "Each box summarises every month in which the zone was observed; zones with shorter records have narrower boxes for that reason alone."
  ) +
  theme_pub(base_size = 9, grid = "x") +
  theme(legend.position = "right")

save_fig(p2, "fig02_distribution_by_zone", DIR_OUT, 9, 12.5, TAG)


# --- Figure 3: monthly climatology heat map ----------------------------------

p3 <- climatology |>
  dplyr::mutate(month_abb = factor(month.abb[month], levels = month.abb),
                zone_name = stats::reorder(zone_name, lst_month_mean,
                                           FUN = mean)) |>
  ggplot(aes(x = month_abb, y = zone_name, fill = lst_month_mean)) +
  geom_tile(colour = "white", linewidth = 0.12) +
  scale_fill_gradientn(colours = PAL_LST, name = "Mean\n(\u00b0C)") +
  scale_x_discrete(expand = c(0, 0)) +
  labs(
    title    = sprintf("Seasonal cycle of surface temperature by zone, %s", AOI_LABEL),
    subtitle = "Mean value of each zone in each calendar month, pooling all years. Horizontal banding indicates that a zone keeps its relative position across seasons.",
    x = NULL, y = NULL,
    caption = "Vertical banding reflects the seasonal cycle shared by the whole city; horizontal banding reflects persistent differences between zones."
  ) +
  theme_pub(base_size = 9, grid = "none") +
  theme(legend.position = "right")

save_fig(p3, "fig03_climatology_heatmap", DIR_OUT, 9, 12.5, TAG)


# --- Figure 4: interannual heat map ------------------------------------------

p4 <- annual |>
  dplyr::mutate(zone_name = stats::reorder(zone_name, lst_year_mean,
                                           FUN = mean)) |>
  ggplot(aes(x = factor(year), y = zone_name, fill = lst_year_mean)) +
  geom_tile(colour = "white", linewidth = 0.12) +
  scale_fill_gradientn(colours = PAL_LST, name = "Mean\n(\u00b0C)",
                       na.value = "grey92") +
  scale_x_discrete(expand = c(0, 0)) +
  labs(
    title    = sprintf("Annual mean surface temperature by zone, %s", AOI_LABEL),
    subtitle = "Blank cells are zone-years with no valid observation. Whole pale columns correspond to years with sparse satellite coverage.",
    x = NULL, y = NULL,
    caption = "Annual means from years with few contributing months are not comparable to well-sampled years; see the coverage calendar in 01_coverage_calendar."
  ) +
  theme_pub(base_size = 9, grid = "none") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))

save_fig(p4, "fig04_annual_heatmap", DIR_OUT, 9.5, 12.5, TAG)


# --- Figure 5: temporal variability by zone ----------------------------------

p5 <- zone_stats |>
  dplyr::arrange(lst_sd_time) |>
  dplyr::mutate(zone_name = factor(zone_name, levels = zone_name)) |>
  ggplot(aes(x = lst_sd_time, y = zone_name, fill = lst_sd_time)) +
  geom_col(width = 0.75) +
  scale_fill_gradientn(colours = PAL_LST, name = "sd (\u00b0C)") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.03))) +
  labs(
    title    = sprintf("Temporal variability of surface temperature by zone, %s", AOI_LABEL),
    subtitle = "Standard deviation of the monthly series within each zone. Larger values mean the zone's apparent temperature swings more between months.",
    x = "Standard deviation across months (\u00b0C)", y = NULL,
    caption = "Part of this variability is the shared seasonal cycle; the ranking is informative, the absolute level is not a pure zone property."
  ) +
  theme_pub(base_size = 9, grid = "x") +
  theme(legend.position = "right")

save_fig(p5, "fig05_variability_by_zone", DIR_OUT, 8.5, 12.5, TAG)


# --- Figure 6: long-term trend by zone ---------------------------------------

trend_df <- zone_stats |> dplyr::filter(!is.na(trend_c_per_year))

if (nrow(trend_df) > 0) {
  p6 <- trend_df |>
    dplyr::arrange(trend_c_per_year) |>
    dplyr::mutate(zone_name = factor(zone_name, levels = zone_name)) |>
    ggplot(aes(x = trend_c_per_year, y = zone_name, fill = trend_c_per_year)) +
    geom_col(width = 0.75) +
    geom_vline(xintercept = 0, colour = "grey25", linewidth = 0.4) +
    scale_fill_gradient2(low = "#2c7bb6", mid = "grey95", high = "#c0392b",
                         midpoint = 0, name = "\u00b0C / year") +
    labs(
      title    = sprintf("Linear trend in surface temperature by zone, %s (%d-%d)",
                         AOI_LABEL, YEAR_START, YEAR_END),
      subtitle = "Slope of a regression on decimal year. Positive values indicate warming over the study period.",
      x = "Trend (\u00b0C per year)", y = NULL,
      caption = sprintf(
        "Zones with fewer than %d monthly observations are omitted. These are descriptive slopes with no significance testing and no adjustment for the uneven seasonal sampling that affects several years.",
        MIN_MONTHS_FOR_TREND)
    ) +
    theme_pub(base_size = 9, grid = "x") +
    theme(legend.position = "right")

  save_fig(p6, "fig06_trend_by_zone", DIR_OUT, 8.5, 12.5, TAG)
}


message(sprintf("%s Done: %s", TAG, Sys.time()))
message(sprintf("%s Outputs in %s\n", TAG, DIR_OUT))
