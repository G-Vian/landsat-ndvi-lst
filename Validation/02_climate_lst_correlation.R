# =============================================================================
# VALIDATION - 02 - SATELLITE LST AGAINST GROUND AIR TEMPERATURE
# =============================================================================
# QUESTION ANSWERED
#   "Does the satellite-derived land surface temperature track the air
#    temperature measured by the ground weather station?"
#
# READ THIS BEFORE INTERPRETING ANY OUTPUT
#   LST and air temperature are NOT the same quantity and are not expected to
#   agree in absolute value. LST is the radiometric temperature of the SURFACE
#   itself - roof tiles, asphalt, grass canopy - retrieved from thermal
#   emission. Air temperature is measured in a shaded screen roughly two metres
#   above ground. On a clear summer day an urban surface can sit 8-15 C above
#   the air around it, so a systematic positive offset is the expected result,
#   not a defect.
#
#   What this script therefore tests is not equality but COVARIATION: do the
#   two series rise and fall together, do they share the same seasonal cycle,
#   and is the relationship stable enough that the satellite record can be
#   treated as a meaningful thermal signal. A high correlation with a large
#   positive offset is exactly the outcome that validates the retrieval; a low
#   correlation would indicate a problem regardless of how close the means are.
#
# WHAT THE SCRIPT DOES
#   1. Reads the monthly LST tables produced by the pipeline.
#   2. Reads the monthly climate spreadsheet and reports exactly which columns
#      it used, so a mismatched column name is caught immediately.
#   3. Joins the two on (year, month) and excludes flagged months.
#   4. Computes Pearson and Spearman correlations and fits LST ~ air
#      temperature for the city and, separately, for every zone.
#   5. Writes CSV tables, an interpreted text report and five figures.
#
# HOW TO RUN
#   singularity exec --bind /home/g.vian:/home/g.vian \
#     /home/public/R_inla/r_inla.sif \
#     Rscript 02_climate_lst_correlation.R
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

suppressPackageStartupMessages(library(readxl))

TAG     <- "[02_climate]"
DIR_OUT <- file.path(VALIDATION_ROOT, "02_climate_lst")
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# CONFIGURATION  ### EDIT HERE ###
# -----------------------------------------------------------------------------
# Monthly climate spreadsheet. The first two columns are assumed to be year and
# month; everything else is matched by the column names below.
PATH_CLIMATE <- file.path(
  "/home/g.vian/Pesquisa_Epidemic/PROJETO_SANTOS",
  "Clima_mensal", "Report_Climate_Santos_2010_2025.xlsx"
)

# Column names inside that spreadsheet, copied EXACTLY as they appear including
# accents, spaces and parentheses. If the script stops with a missing-column
# error it prints every name it did find, so the fix is to copy the right one
# from that list into the corresponding line here.
COL_TEMP_MEAN <- "Mean Temperature (\u00b0C)"
COL_TEMP_MIN  <- "Min Temperature (\u00b0C)"
COL_TEMP_MAX  <- "Max Temperature (\u00b0C)"
# -----------------------------------------------------------------------------

message(sprintf("\n%s Start: %s", TAG, Sys.time()))
message(sprintf("%s Output: %s", TAG, DIR_OUT))


# -----------------------------------------------------------------------------
# 1. LOAD SATELLITE DATA
# -----------------------------------------------------------------------------
pipeline_dir <- resolve_pipeline_dir()
tables_dir   <- resolve_tables_dir(pipeline_dir)

city <- load_city_monthly(tables_dir) |>
  dplyr::select(dplyr::any_of(c("year", "month", "lst_mean", "lst_min",
                                "lst_max", "ndvi_mean", "n_scenes")))

zones <- tryCatch(load_zone_monthly(tables_dir), error = function(e) NULL)


# -----------------------------------------------------------------------------
# 2. LOAD CLIMATE DATA
# -----------------------------------------------------------------------------
if (!file.exists(PATH_CLIMATE))
  stop("Climate spreadsheet not found:\n  ", PATH_CLIMATE,
       "\nEdit PATH_CLIMATE at the top of this script.")

climate_raw <- readxl::read_excel(PATH_CLIMATE)
names(climate_raw)[1:2] <- c("year", "month")

# Transparency log. Silently picking the wrong temperature column would produce
# a plausible-looking but meaningless correlation, so what was read is printed
# before anything is computed.
message(sprintf("%s Climate file: %s", TAG, basename(PATH_CLIMATE)))
message(sprintf("%s   rows: %d | period: %s-%s", TAG, nrow(climate_raw),
                min(climate_raw$year, na.rm = TRUE),
                max(climate_raw$year, na.rm = TRUE)))
message(sprintf("%s   columns present:", TAG))
for (cn in names(climate_raw)) message(sprintf("%s     - %s", TAG, cn))

# Coerce every column to numeric, converting decimal commas to points. A single
# text cell anywhere in the sheet would otherwise make the whole column
# character, and the correlation would fail with an opaque error.
climate_raw <- climate_raw |>
  dplyr::mutate(dplyr::across(
    dplyr::everything(),
    ~ suppressWarnings(as.numeric(gsub(",", ".", as.character(.))))))

needed  <- c(COL_TEMP_MEAN, COL_TEMP_MIN, COL_TEMP_MAX)
missing <- setdiff(needed, names(climate_raw))
if (length(missing) > 0)
  stop("Column(s) not found in the climate spreadsheet:\n  ",
       paste(missing, collapse = "\n  "),
       "\n\nColumns that ARE present:\n  ",
       paste(names(climate_raw), collapse = "\n  "),
       "\n\nCopy the correct names into COL_TEMP_* at the top of this script.")

climate <- climate_raw |>
  dplyr::transmute(
    year, month,
    temp_mean = .data[[COL_TEMP_MEAN]],
    temp_min  = .data[[COL_TEMP_MIN]],
    temp_max  = .data[[COL_TEMP_MAX]]
  )

describe <- function(x, label) {
  n_na <- sum(is.na(x))
  if (n_na == length(x)) {
    message(sprintf("%s     %-10s : ALL MISSING (%d values)", TAG, label, n_na))
  } else {
    message(sprintf("%s     %-10s : mean %.2f | range %.2f to %.2f | %d missing",
                    TAG, label, mean(x, na.rm = TRUE), min(x, na.rm = TRUE),
                    max(x, na.rm = TRUE), n_na))
  }
}
message(sprintf("%s   summary of the three columns used:", TAG))
describe(climate$temp_mean, "temp_mean")
describe(climate$temp_min,  "temp_min")
describe(climate$temp_max,  "temp_max")


# -----------------------------------------------------------------------------
# 3. JOIN AND FLAG
# -----------------------------------------------------------------------------
# An inner join is deliberate: a month is only usable when BOTH the satellite
# and the station observed it. Flagged months are kept in the data frame so
# they can be drawn on the figures, but they are removed from every statistic.

df <- city |>
  dplyr::inner_join(climate, by = c("year", "month")) |>
  add_time_columns()

df_clean <- df |>
  dplyr::filter(!low_quality, !is.na(lst_mean), !is.na(temp_mean))

message(sprintf("%s %d months with both satellite and station data (after flags).",
                TAG, nrow(df_clean)))

if (nrow(df_clean) < 10)
  stop("Fewer than 10 usable months. Check that the climate spreadsheet and ",
       "the satellite tables cover the same period.")


# -----------------------------------------------------------------------------
# 4. CORRELATION AND REGRESSION, CITY LEVEL
# -----------------------------------------------------------------------------
# Both a parametric and a rank correlation are reported. Pearson measures
# linear association and is the quantity the regression slope refers to;
# Spearman is insensitive to outliers and to any monotone nonlinearity, so a
# large gap between the two would itself be diagnostic.

corr_pair <- function(x, y) {
  ok <- stats::complete.cases(x, y)
  c(pearson  = stats::cor(x[ok], y[ok], method = "pearson"),
    spearman = stats::cor(x[ok], y[ok], method = "spearman"),
    n        = sum(ok))
}

c_mean <- corr_pair(df_clean$temp_mean, df_clean$lst_mean)
c_min  <- corr_pair(df_clean$temp_min,  df_clean$lst_mean)
c_max  <- corr_pair(df_clean$temp_max,  df_clean$lst_mean)

fit_mean <- stats::lm(lst_mean ~ temp_mean, data = df_clean)
fit_min  <- stats::lm(lst_mean ~ temp_min,  data = df_clean)
fit_max  <- stats::lm(lst_mean ~ temp_max,  data = df_clean)

corr_table <- dplyr::tibble(
  climate_variable = c("Mean air temperature", "Minimum air temperature",
                       "Maximum air temperature"),
  pearson   = round(c(c_mean["pearson"],  c_min["pearson"],  c_max["pearson"]), 3),
  spearman  = round(c(c_mean["spearman"], c_min["spearman"], c_max["spearman"]), 3),
  n_months  = c(c_mean["n"], c_min["n"], c_max["n"]),
  intercept = round(c(stats::coef(fit_mean)[1], stats::coef(fit_min)[1],
                      stats::coef(fit_max)[1]), 3),
  slope     = round(c(stats::coef(fit_mean)[2], stats::coef(fit_min)[2],
                      stats::coef(fit_max)[2]), 3),
  r_squared = round(c(summary(fit_mean)$r.squared, summary(fit_min)$r.squared,
                      summary(fit_max)$r.squared), 3)
)
readr::write_csv(corr_table, file.path(DIR_OUT, "correlation_city.csv"))


# -----------------------------------------------------------------------------
# 5. REGRESSION PER ZONE
# -----------------------------------------------------------------------------
# Fitting the same relationship separately in every zone shows whether the
# agreement is uniform across the city or concentrated in particular areas.
# A zone whose surface is largely vegetated tracks the air closely; a densely
# built one has a much steeper response, which is the urban heat island
# expressing itself.

zone_fit <- NULL
if (!is.null(zones)) {

  zone_join <- zones |>
    dplyr::select(zone_name, year, month, lst_mean) |>
    dplyr::inner_join(climate, by = c("year", "month")) |>
    add_time_columns() |>
    dplyr::filter(!low_quality, !is.na(lst_mean), !is.na(temp_mean))

  # A slope fitted to a handful of months is dominated by which months happen
  # to be present, so zones with a thin record are dropped rather than reported
  # with a spuriously precise coefficient.
  MIN_MONTHS_PER_ZONE <- 12

  zone_fit <- zone_join |>
    dplyr::group_by(zone_name) |>
    dplyr::filter(dplyr::n() >= MIN_MONTHS_PER_ZONE) |>
    dplyr::summarise(
      n_months  = dplyr::n(),
      pearson   = stats::cor(temp_mean, lst_mean, method = "pearson"),
      spearman  = stats::cor(temp_mean, lst_mean, method = "spearman"),
      intercept = unname(stats::coef(stats::lm(lst_mean ~ temp_mean))[1]),
      slope     = unname(stats::coef(stats::lm(lst_mean ~ temp_mean))[2]),
      lst_overall_mean = mean(lst_mean, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~round(., 3))) |>
    dplyr::arrange(dplyr::desc(pearson))

  readr::write_csv(zone_fit, file.path(DIR_OUT, "regression_by_zone.csv"))
  message(sprintf("%s %d zones with at least %d months.",
                  TAG, nrow(zone_fit), MIN_MONTHS_PER_ZONE))
}


# -----------------------------------------------------------------------------
# 6. TEXT REPORT
# -----------------------------------------------------------------------------
b0 <- stats::coef(fit_mean)[1]
b1 <- stats::coef(fit_mean)[2]
r2 <- summary(fit_mean)$r.squared

report <- file.path(DIR_OUT, "climate_lst_report.txt")
open_report(report)

cat(strrep("=", 78), "\n")
cat("  SATELLITE LST vs GROUND AIR TEMPERATURE -", AOI_LABEL, "\n")
cat(sprintf("  Generated: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(strrep("=", 78), "\n\n")

cat("--- DATA SOURCES ---\n")
cat(sprintf("  Climate file : %s\n", PATH_CLIMATE))
cat(sprintf("  Mean column  : '%s'\n", COL_TEMP_MEAN))
cat(sprintf("  Min column   : '%s'\n", COL_TEMP_MIN))
cat(sprintf("  Max column   : '%s'\n", COL_TEMP_MAX))
cat(sprintf("  Satellite    : %s\n", pipeline_dir))
cat(sprintf("  Months joined: %d (after removing flagged months)\n\n",
            nrow(df_clean)))

cat("--- WHAT IS AND IS NOT BEING TESTED ---\n")
cat("LST is the radiometric temperature of the surface itself; air temperature\n")
cat("is measured in a shaded screen about two metres above ground. On a clear\n")
cat("summer day an urban surface commonly sits 8-15 C above the surrounding\n")
cat("air, so a positive offset is the expected result and not a defect.\n")
cat("The test here is whether the two series COVARY, not whether they match.\n\n")

cat("--- CITY-LEVEL CORRELATION AND REGRESSION ---\n")
print(as.data.frame(corr_table), row.names = FALSE)
cat("\n")

cat("--- INTERPRETATION ---\n")
cat(sprintf("  LST = %.2f + %.3f x (mean air temperature)     R2 = %.3f\n\n",
            b0, b1, r2))
if (b1 > 0) {
  cat(sprintf("  Each 1 C rise in air temperature corresponds to a %.2f C rise\n", b1))
  cat("  in the city-wide surface temperature.\n")
  if (b1 > 1.2) {
    cat("  A slope above one means the surface amplifies the atmospheric signal,\n")
    cat("  which is the characteristic behaviour of a built-up surface with low\n")
    cat("  evaporative cooling.\n")
  } else if (b1 < 0.8) {
    cat("  A slope below one means the surface damps the atmospheric signal,\n")
    cat("  consistent with substantial vegetation or water in the footprint.\n")
  }
} else {
  cat("  NEGATIVE slope. This is not physically expected and points to a data\n")
  cat("  problem: check that the climate columns were read correctly and that\n")
  cat("  the year/month join is aligned.\n")
}
cat("\n")

r_abs <- abs(c_mean["pearson"])
if (r_abs > 0.7) {
  cat("  Correlation is STRONG (|r| > 0.7). The satellite series follows the\n")
  cat("  same seasonal regime as the meteorological record, which is the\n")
  cat("  expected outcome and supports the retrieval.\n")
} else if (r_abs > 0.4) {
  cat("  Correlation is MODERATE (0.4 < |r| < 0.7). The seasonal cycle is\n")
  cat("  shared but with substantial scatter; inspect monthly coverage and the\n")
  cat("  flagged months before drawing conclusions.\n")
} else {
  cat("  Correlation is WEAK. Investigate: likely causes are a misaligned\n")
  cat("  year/month join, the wrong climate column, or unfiltered artefacts.\n")
}
cat("\n")

diffs <- df_clean$lst_mean - df_clean$temp_mean
cat("--- SURFACE-TO-AIR OFFSET ---\n")
cat(sprintf("  mean %.1f C | median %.1f C | sd %.1f C | range %.1f to %.1f C\n",
            mean(diffs), stats::median(diffs), stats::sd(diffs),
            min(diffs), max(diffs)))
cat("  A consistently positive offset of roughly this size is what the urban\n")
cat("  heat island literature reports for daytime overpasses.\n\n")

if (!is.null(zone_fit)) {
  cat("--- HETEROGENEITY ACROSS ZONES ---\n")
  cat(sprintf("  %d zones analysed.\n", nrow(zone_fit)))
  cat(sprintf("  Pearson r : median %.3f | range %.3f to %.3f\n",
              stats::median(zone_fit$pearson, na.rm = TRUE),
              min(zone_fit$pearson, na.rm = TRUE),
              max(zone_fit$pearson, na.rm = TRUE)))
  cat(sprintf("  Slope     : median %.3f | range %.3f to %.3f\n\n",
              stats::median(zone_fit$slope, na.rm = TRUE),
              min(zone_fit$slope, na.rm = TRUE),
              max(zone_fit$slope, na.rm = TRUE)))
  cat("  Zones with the strongest coupling to air temperature:\n")
  print(as.data.frame(utils::head(
    zone_fit[, c("zone_name", "n_months", "pearson", "slope")], 5)),
    row.names = FALSE)
  cat("\n  Zones with the weakest coupling:\n")
  print(as.data.frame(utils::tail(
    zone_fit[, c("zone_name", "n_months", "pearson", "slope")], 5)),
    row.names = FALSE)
  cat("\n")
}

cat("--- EXCLUSIONS ---\n")
for (i in seq_len(nrow(LOW_QUALITY_MONTHS)))
  cat(sprintf("  %d-%02d excluded: %s\n", LOW_QUALITY_MONTHS$year[i],
              LOW_QUALITY_MONTHS$month[i], LOW_QUALITY_MONTHS$reason[i]))
cat("\n")
cat(strrep("=", 78), "\n")

close_report()
message(sprintf("%s Report: %s", TAG, basename(report)))


# -----------------------------------------------------------------------------
# 7. FIGURES
# -----------------------------------------------------------------------------

# --- Figure 1: scatter with regression, coloured by season -------------------
# Colouring by season shows whether the relationship is driven purely by the
# annual cycle (points would form a single elongated cloud ordered by season)
# or whether it also holds within seasons (spread within each colour).

p1 <- ggplot(df_clean, aes(x = temp_mean, y = lst_mean)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "grey25", fill = "grey80", linewidth = 0.7) +
  geom_point(aes(fill = season), shape = 21, size = 3,
             colour = "white", stroke = 0.5, alpha = 0.95) +
  # Flagged months are drawn as open crosses so the reader can see they exist
  # and can see that they sit off the trend, without them entering the fit.
  geom_point(data = dplyr::filter(df, low_quality),
             shape = 4, size = 3.4, stroke = 1.1, colour = "#c0392b") +
  scale_fill_manual(values = PAL_SEASON, name = "Season",
                    limits = c("Summer", "Autumn", "Winter", "Spring")) +
  labs(
    title    = sprintf("Surface temperature against air temperature, %s", AOI_LABEL),
    subtitle = sprintf(
      "Pearson r = %.2f | LST = %.1f + %.2f x air temperature | R2 = %.2f | n = %d months",
      c_mean["pearson"], b0, b1, r2, nrow(df_clean)),
    x = "Mean air temperature (\u00b0C), weather station",
    y = "Mean land surface temperature (\u00b0C), Landsat",
    caption = paste0(
      "Line is an ordinary least squares fit with a 95% confidence band. ",
      "Red crosses mark flagged months, excluded from the fit.\n",
      "Surface temperature exceeding air temperature is expected for a ",
      "mid-morning overpass over an urban surface.")
  ) +
  theme_pub() +
  theme(legend.position = "right")

save_fig(p1, "fig01_scatter_lst_vs_air", DIR_OUT, 8.5, 6, TAG)


# --- Figure 2: overlaid time series ------------------------------------------
# The point of this panel is the PARALLELISM of the curves, not their overlap.
# A constant vertical gap between them is the surface-to-air offset.

series <- df |>
  dplyr::select(date, lst_mean, temp_mean, temp_max, low_quality) |>
  tidyr::pivot_longer(c(lst_mean, temp_mean, temp_max),
                      names_to = "variable", values_to = "value") |>
  dplyr::mutate(variable = factor(dplyr::recode(
    variable,
    lst_mean  = "Surface temperature (Landsat)",
    temp_mean = "Mean air temperature (station)",
    temp_max  = "Maximum air temperature (station)"),
    levels = c("Surface temperature (Landsat)",
               "Maximum air temperature (station)",
               "Mean air temperature (station)")))

p2 <- ggplot(series, aes(x = date, y = value, colour = variable)) +
  geom_line(linewidth = 0.5, alpha = 0.9) +
  geom_point(size = 1.1, alpha = 0.9) +
  scale_colour_manual(
    values = c("Surface temperature (Landsat)"     = "#c0392b",
               "Maximum air temperature (station)" = "#e08214",
               "Mean air temperature (station)"    = "#2c7bb6"),
    name = NULL) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y",
               expand = expansion(mult = 0.01)) +
  guides(colour = guide_legend(nrow = 1)) +
  labs(
    title    = sprintf("Monthly time series of surface and air temperature, %s", AOI_LABEL),
    subtitle = "The diagnostic feature is that the curves move in parallel; the vertical gap between them is the surface-to-air offset.",
    x = NULL, y = "Temperature (\u00b0C)",
    caption = "Gaps in the satellite curve are months without a cloud-free acquisition."
  ) +
  theme_pub(grid = "y") +
  theme(legend.position = "bottom")

save_fig(p2, "fig02_timeseries_overlay", DIR_OUT, 11, 5.5, TAG)


# --- Figure 3: monthly climatology ------------------------------------------
# Side-by-side boxplots per calendar month test whether the two records agree
# on the SHAPE of the seasonal cycle, independently of their level.

box_df <- df_clean |>
  dplyr::select(month_abb, lst_mean, temp_mean) |>
  tidyr::pivot_longer(c(lst_mean, temp_mean),
                      names_to = "variable", values_to = "value") |>
  dplyr::mutate(variable = dplyr::recode(
    variable,
    lst_mean  = "Surface (Landsat)",
    temp_mean = "Air (station)"))

p3 <- ggplot(box_df, aes(x = month_abb, y = value, fill = variable)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.62,
               alpha = 0.9, outlier.size = 0.9, linewidth = 0.35) +
  scale_fill_manual(values = c("Surface (Landsat)" = "#c0392b",
                               "Air (station)"     = "#2c7bb6"),
                    name = NULL) +
  labs(
    title    = sprintf("Seasonal cycle of surface and air temperature, %s", AOI_LABEL),
    subtitle = "Both records peak in the austral summer and trough in winter; the offset between them widens in the warm months.",
    x = NULL, y = "Temperature (\u00b0C)",
    caption = sprintf("Each box aggregates every occurrence of that calendar month between %d and %d.",
                      YEAR_START, YEAR_END)
  ) +
  theme_pub(grid = "y") +
  theme(legend.position = "bottom")

save_fig(p3, "fig03_monthly_climatology", DIR_OUT, 10, 5.5, TAG)


# --- Figure 4: distribution of the surface-to-air offset ---------------------

diff_df <- df_clean |> dplyr::mutate(offset = lst_mean - temp_mean)

p4 <- ggplot(diff_df, aes(x = offset)) +
  geom_histogram(bins = 28, fill = "#e08214", colour = "white",
                 linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = "grey25", linewidth = 0.5) +
  geom_vline(xintercept = mean(diff_df$offset), colour = "#c0392b",
             linewidth = 0.9, linetype = "dashed") +
  annotate("text", x = mean(diff_df$offset), y = Inf,
           label = sprintf("  mean %.1f \u00b0C", mean(diff_df$offset)),
           hjust = 0, vjust = 1.8, size = 3.2, colour = "#c0392b") +
  labs(
    title    = sprintf("Surface-to-air temperature offset, %s", AOI_LABEL),
    subtitle = sprintf("mean %.1f \u00b0C | median %.1f \u00b0C | sd %.1f \u00b0C | n = %d months",
                       mean(diff_df$offset), stats::median(diff_df$offset),
                       stats::sd(diff_df$offset), nrow(diff_df)),
    x = "Surface minus air temperature (\u00b0C)", y = "Number of months",
    caption = paste0(
      "Values to the right of the solid line mean the surface is warmer than ",
      "the air, which is the expected daytime condition.\n",
      "A distribution centred near zero or extending below it would suggest a ",
      "retrieval or calibration problem.")
  ) +
  theme_pub(grid = "y")

save_fig(p4, "fig04_surface_air_offset", DIR_OUT, 8, 5.5, TAG)


# --- Figure 5: per-zone slope and correlation --------------------------------
# Two horizontal bar charts side by side would be wide and hard to read at
# 60 zones, so they are produced as separate tall figures.

if (!is.null(zone_fit) && nrow(zone_fit) > 0) {

  p5 <- zone_fit |>
    dplyr::arrange(slope) |>
    dplyr::mutate(zone_name = factor(zone_name, levels = zone_name)) |>
    ggplot(aes(x = slope, y = zone_name, fill = slope)) +
    geom_col(width = 0.75) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey30",
               linewidth = 0.4) +
    scale_fill_gradientn(colours = PAL_LST, name = "Slope\n(\u00b0C / \u00b0C)") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(
      title    = "Sensitivity of surface temperature to air temperature, by zone",
      subtitle = "Regression slope of LST on air temperature. The dashed line at one separates zones that amplify the atmospheric signal from those that damp it.",
      x = "Slope (\u00b0C surface per \u00b0C air)", y = NULL,
      caption = sprintf("Zones with fewer than 12 monthly observations are omitted. n = %d zones.",
                        nrow(zone_fit))
    ) +
    theme_pub(base_size = 9, grid = "x") +
    theme(legend.position = "right")

  save_fig(p5, "fig05_slope_by_zone", DIR_OUT, 8.5, 12, TAG)

  p6 <- zone_fit |>
    dplyr::arrange(pearson) |>
    dplyr::mutate(zone_name = factor(zone_name, levels = zone_name)) |>
    ggplot(aes(x = pearson, y = zone_name, fill = pearson)) +
    geom_col(width = 0.75) +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "grey30",
               linewidth = 0.4) +
    scale_fill_gradientn(colours = PAL_LST, limits = c(-1, 1),
                         name = "Pearson r") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(
      title    = "Agreement between surface and air temperature, by zone",
      subtitle = "Pearson correlation of the monthly series. The dashed line marks the conventional threshold for a strong association.",
      x = "Pearson correlation", y = NULL,
      caption = "A zone falling well below the others usually has a short or patchy record; cross-check n_months in regression_by_zone.csv."
    ) +
    theme_pub(base_size = 9, grid = "x") +
    theme(legend.position = "right")

  save_fig(p6, "fig06_correlation_by_zone", DIR_OUT, 8.5, 12, TAG)
}


message(sprintf("%s Done: %s", TAG, Sys.time()))
message(sprintf("%s Outputs in %s\n", TAG, DIR_OUT))
