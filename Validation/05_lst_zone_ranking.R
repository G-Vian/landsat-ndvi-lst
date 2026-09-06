# =============================================================================
# VALIDATION - 05 - PERSISTENCE OF THE SURFACE TEMPERATURE RANKING
# =============================================================================
# QUESTION ANSWERED
#   "Are the same zones always the warmest, or the coolest, regardless of month
#    and year? In other words, is the spatial pattern of surface temperature a
#    stable property of the city, or does it reshuffle over time?"
#
# WHY THIS IS THE MOST IMPORTANT OF THE FIVE SCRIPTS FOR A SPATIAL MODEL
#   A zone-level model attributes an exposure to each zone. That is only
#   defensible if the ordering of zones is stable: if the warmest zone changes
#   every month, then a zone's long-run mean temperature is not describing a
#   property of the place, it is averaging noise. Conversely, if a handful of
#   zones sit at the top of the ranking in almost every month of a sixteen-year
#   record, the spatial pattern is structural and the exposure is meaningful.
#
#   Note that this is a question about RANK, not about level. The absolute
#   temperature of every zone rises in summer and falls in winter; that shared
#   seasonal movement is removed by ranking within each month, leaving only the
#   relative positions.
#
# METHOD
#   Within each (year, month) the zones are ranked from coolest (rank 1) to
#   warmest (rank N). Each month is effectively an independent vote on the
#   spatial ordering. Per zone the script then reports:
#     mean_rank        average position across the record
#     sd_rank          how much the position moves; low means persistent
#     pct_top_warm     share of months spent among the N warmest
#     pct_top_cool     share of months spent among the N coolest
#     delta_vs_city    mean difference from the city-wide average of the month
#     r_vs_city        correlation of the zone series with the city series
#
# HOW TO RUN
#   singularity exec --bind $HOME:$HOME \
#     /path/to/container.sif \
#     Rscript 05_lst_zone_ranking.R
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

TAG     <- "[05_ranking]"
DIR_OUT <- file.path(VALIDATION_ROOT, "05_lst_ranking")
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

# How many zones count as "the top". With about 60 zones, 10 is roughly the top
# and bottom sixth, which is a natural definition of an extreme group.
N_TOP <- 10

message(sprintf("\n%s Start: %s", TAG, Sys.time()))
message(sprintf("%s Output: %s", TAG, DIR_OUT))


# -----------------------------------------------------------------------------
# 1. LOAD
# -----------------------------------------------------------------------------
pipeline_dir <- resolve_pipeline_dir()
tables_dir   <- resolve_tables_dir(pipeline_dir)

df <- load_zone_monthly(tables_dir) |>
  dplyr::select(dplyr::any_of(c("zone_id", "zone_name", "year", "month",
                                "lst_mean"))) |>
  add_time_columns() |>
  dplyr::filter(!is.na(lst_mean), !low_quality)

# Resolve the study period. YEAR_START / YEAR_END may be left NULL in
# 00_common.R, in which case the span is taken from the data itself, so the
# report always describes the run that actually happened.
resolve_period(df)

n_zones  <- dplyr::n_distinct(df$zone_name)
n_months <- dplyr::n_distinct(df$date)
message(sprintf("%s %d zones, %d months with data.", TAG, n_zones, n_months))


# -----------------------------------------------------------------------------
# 2. RANK WITHIN EACH MONTH
# -----------------------------------------------------------------------------
# Ranking inside the month is what removes the shared seasonal cycle. Ties get
# the average rank, which matters when several zones return identical rounded
# values.
#
# rank_norm rescales the rank to [0, 1] so that months in which fewer zones
# were observed remain comparable with fully observed months. Without it, a
# month with 30 observed zones would put its warmest zone at rank 30 while a
# complete month puts its warmest at 60, and the two would not be comparable.

ranked <- df |>
  dplyr::group_by(year, month) |>
  dplyr::mutate(
    n_zones_month = dplyr::n(),
    rank      = rank(lst_mean, ties.method = "average"),
    rank_norm = (rank - 1) / (n_zones_month - 1),
    group     = dplyr::case_when(
      rank <= N_TOP                    ~ "coolest",
      rank >  n_zones_month - N_TOP    ~ "warmest",
      TRUE                             ~ "middle"
    )
  ) |>
  dplyr::ungroup()

readr::write_csv(
  dplyr::select(ranked, zone_name, year, month, date, lst_mean,
                rank, rank_norm, group, n_zones_month),
  file.path(DIR_OUT, "ranking_by_month.csv"))


# -----------------------------------------------------------------------------
# 3. DIFFERENCE FROM THE CITY MEAN
# -----------------------------------------------------------------------------
# A zone that sits above the city average in nearly every month is a stable
# heat island; one that oscillates around it is not. This is the same anomaly
# construction used to build the model covariate, so the numbers here describe
# the covariate directly.

city_by_month <- df |>
  dplyr::group_by(year, month) |>
  dplyr::summarise(lst_city = mean(lst_mean, na.rm = TRUE), .groups = "drop")

with_city <- df |>
  dplyr::left_join(city_by_month, by = c("year", "month")) |>
  dplyr::mutate(delta = lst_mean - lst_city)


# -----------------------------------------------------------------------------
# 4. PERSISTENCE STATISTICS PER ZONE
# -----------------------------------------------------------------------------
persistence <- ranked |>
  dplyr::left_join(dplyr::select(with_city, zone_name, year, month, delta),
                   by = c("zone_name", "year", "month")) |>
  dplyr::group_by(zone_id, zone_name) |>
  dplyr::summarise(
    n_months_obs   = dplyr::n(),
    lst_mean_all   = mean(lst_mean, na.rm = TRUE),
    mean_rank      = mean(rank, na.rm = TRUE),
    sd_rank        = stats::sd(rank, na.rm = TRUE),
    min_rank       = min(rank, na.rm = TRUE),
    max_rank       = max(rank, na.rm = TRUE),
    mean_rank_norm = mean(rank_norm, na.rm = TRUE),
    pct_top_warm   = round(100 * mean(group == "warmest", na.rm = TRUE), 1),
    pct_top_cool   = round(100 * mean(group == "coolest", na.rm = TRUE), 1),
    delta_vs_city  = mean(delta, na.rm = TRUE),
    pct_above_city = round(100 * mean(delta > 0, na.rm = TRUE), 1),
    .groups = "drop"
  )

# Correlation of each zone series with the city series. A value near one means
# the zone simply follows the city up and down; a lower value means the zone
# has its own dynamics on top of the shared seasonal cycle.
corr_city <- with_city |>
  dplyr::group_by(zone_name) |>
  dplyr::summarise(
    n_pairs   = sum(!is.na(lst_mean) & !is.na(lst_city)),
    r_vs_city = if (sum(!is.na(lst_mean) & !is.na(lst_city)) >= 3)
      stats::cor(lst_mean, lst_city, use = "complete.obs") else NA_real_,
    .groups = "drop"
  )

persistence <- persistence |>
  dplyr::left_join(dplyr::select(corr_city, zone_name, r_vs_city),
                   by = "zone_name") |>
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~round(., 3))) |>
  dplyr::arrange(dplyr::desc(mean_rank))

readr::write_csv(persistence,
                 file.path(DIR_OUT, "persistence_statistics.csv"))

# Spatial spread within each month: how far apart are the warmest and coolest
# zones on any given date.
amplitude <- ranked |>
  dplyr::group_by(year, month, date) |>
  dplyr::summarise(
    amplitude = max(lst_mean, na.rm = TRUE) - min(lst_mean, na.rm = TRUE),
    n_zones_month = dplyr::first(n_zones_month),
    .groups = "drop"
  )
readr::write_csv(amplitude, file.path(DIR_OUT, "spatial_amplitude.csv"))


# -----------------------------------------------------------------------------
# 5. TEXT REPORT
# -----------------------------------------------------------------------------
always_warm <- persistence |> dplyr::filter(pct_top_warm >= 80) |>
  dplyr::arrange(dplyr::desc(pct_top_warm))
always_cool <- persistence |> dplyr::filter(pct_top_cool >= 80) |>
  dplyr::arrange(dplyr::desc(pct_top_cool))
most_stable <- utils::head(dplyr::arrange(persistence, sd_rank), 10)
most_mobile <- utils::head(dplyr::arrange(persistence,
                                          dplyr::desc(sd_rank)), 10)

report <- file.path(DIR_OUT, "ranking_persistence_report.txt")
open_report(report)

cat(strrep("=", 78), "\n")
cat("  PERSISTENCE OF THE SURFACE TEMPERATURE RANKING -", AOI_LABEL, "\n")
cat(sprintf("  Zones: %d | Months: %d | Top group size: %d\n",
            n_zones, n_months, N_TOP))
cat(sprintf("  Generated: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(strrep("=", 78), "\n\n")

cat("--- HOW TO READ THIS REPORT ---\n")
cat("Within each month the zones are ordered by mean surface temperature: the\n")
cat("coolest receives rank 1 and the warmest rank N. Repeating this for every\n")
cat("month of the record and asking who occupies which position separates a\n")
cat("structural spatial pattern from month-to-month noise.\n")
cat("  low  sd_rank -> the zone holds its position -> persistent\n")
cat("  high sd_rank -> the zone moves around the ranking -> volatile\n")
cat(sprintf("  pct_top_warm -> share of months spent among the %d warmest\n\n",
            N_TOP))

cat("--- SPATIAL SPREAD WITHIN A MONTH ---\n")
cat("  Difference between the warmest and coolest zone on the same date:\n")
cat(sprintf("    mean   : %.2f C\n",   mean(amplitude$amplitude)))
cat(sprintf("    median : %.2f C\n",   stats::median(amplitude$amplitude)))
cat(sprintf("    minimum: %.2f C  (%d-%02d)\n", min(amplitude$amplitude),
            amplitude$year[which.min(amplitude$amplitude)],
            amplitude$month[which.min(amplitude$amplitude)]))
cat(sprintf("    maximum: %.2f C  (%d-%02d)\n", max(amplitude$amplitude),
            amplitude$year[which.max(amplitude$amplitude)],
            amplitude$month[which.max(amplitude$amplitude)]))
cat("\n  This is the contrast available to a spatial model on a single date.\n\n")

cat(sprintf("--- ZONES AMONG THE %d WARMEST IN AT LEAST 80%% OF MONTHS ---\n",
            N_TOP))
if (nrow(always_warm) > 0) {
  print(as.data.frame(always_warm[, c("zone_name", "pct_top_warm", "mean_rank",
                                      "sd_rank", "lst_mean_all",
                                      "n_months_obs")]), row.names = FALSE)
} else {
  cat(sprintf("  None. No zone stays in the top %d in 80%% of months.\n", N_TOP))
}

cat(sprintf("\n--- ZONES AMONG THE %d COOLEST IN AT LEAST 80%% OF MONTHS ---\n",
            N_TOP))
if (nrow(always_cool) > 0) {
  print(as.data.frame(always_cool[, c("zone_name", "pct_top_cool", "mean_rank",
                                      "sd_rank", "lst_mean_all",
                                      "n_months_obs")]), row.names = FALSE)
} else {
  cat(sprintf("  None. No zone stays in the bottom %d in 80%% of months.\n",
              N_TOP))
}
cat("\n")

cat("--- TEN ZONES WITH THE MOST STABLE POSITION ---\n")
print(as.data.frame(most_stable[, c("zone_name", "mean_rank", "sd_rank",
                                    "lst_mean_all", "n_months_obs")]),
      row.names = FALSE)
cat("\n--- TEN ZONES WITH THE MOST VOLATILE POSITION ---\n")
print(as.data.frame(most_mobile[, c("zone_name", "mean_rank", "sd_rank",
                                    "lst_mean_all", "n_months_obs")]),
      row.names = FALSE)
cat("\n")

cat("--- COUPLING OF EACH ZONE TO THE CITY AVERAGE ---\n")
cat("  Correlation between a zone's monthly series and the city-wide series.\n")
cat("  Values near one mean the zone simply follows the city up and down.\n\n")
cat(sprintf("  Mean    : %.3f\n", mean(persistence$r_vs_city, na.rm = TRUE)))
cat(sprintf("  Minimum : %.3f (%s)\n",
            min(persistence$r_vs_city, na.rm = TRUE),
            persistence$zone_name[which.min(persistence$r_vs_city)]))
cat(sprintf("  Maximum : %.3f (%s)\n\n",
            max(persistence$r_vs_city, na.rm = TRUE),
            persistence$zone_name[which.max(persistence$r_vs_city)]))

cat("--- CONCLUSION ---\n")
n_warm_half <- sum(persistence$pct_top_warm >= 50, na.rm = TRUE)
n_cool_half <- sum(persistence$pct_top_cool >= 50, na.rm = TRUE)
cat(sprintf("  %d zones sit among the %d warmest in more than half of the\n",
            n_warm_half, N_TOP))
cat("    months observed.\n")
cat(sprintf("  %d zones sit among the %d coolest in more than half of the\n",
            n_cool_half, N_TOP))
cat("    months observed.\n\n")

if (n_warm_half >= 5 && n_cool_half >= 5) {
  cat("  The spatial variation of surface temperature is STRUCTURAL. Some\n")
  cat("  zones are systematically warmer and others systematically cooler,\n")
  cat("  independently of month and year. This is the expected signature of an\n")
  cat("  urban heat island: the absolute temperature follows the season, but\n")
  cat("  the spatial hierarchy underneath it stays fixed.\n\n")
  cat("  For a zone-level statistical model this is the favourable result. It\n")
  cat("  means a zone's long-run temperature describes a property of the place\n")
  cat("  rather than an average of noise, and that a spatial exposure term is\n")
  cat("  justified.\n")
} else {
  cat("  The ranking between zones is NOT stable: few zones hold a consistent\n")
  cat("  position. Treating a zone's long-run mean temperature as a fixed\n")
  cat("  spatial exposure would then be difficult to defend, and a formulation\n")
  cat("  with an explicit space-time interaction should be preferred.\n")
}
cat("\n")
cat(strrep("=", 78), "\n")

close_report()
message(sprintf("%s Report: %s", TAG, basename(report)))


# -----------------------------------------------------------------------------
# 6. FIGURES
# -----------------------------------------------------------------------------

# --- Figure 1: rank heat map through time ------------------------------------
# This is the central figure of the script. Horizontal bands of uniform colour
# running the full width mean a zone kept its position for sixteen years;
# mottled rows mean it moved. Because the colour is the normalised rank rather
# than the temperature, the seasonal cycle has been removed and what remains is
# purely the spatial ordering.

p1 <- ranked |>
  dplyr::mutate(zone_name = stats::reorder(zone_name, rank, FUN = mean)) |>
  ggplot(aes(x = date, y = zone_name, fill = rank_norm)) +
  geom_tile() +
  scale_fill_gradientn(colours = PAL_LST,
                       name = "Rank within\nthe month",
                       breaks = c(0, 0.5, 1),
                       labels = c("Coolest", "Middle", "Warmest")) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y",
               expand = c(0, 0)) +
  labs(
    title    = sprintf("Persistence of the surface temperature ranking, %s", AOI_LABEL),
    subtitle = "Each column is one month; colour is the zone's position within that month, so the shared seasonal cycle is removed.",
    x = NULL, y = NULL,
    caption = paste0(
      "Continuous horizontal bands indicate zones that hold their relative position across the whole record; ",
      "mottled rows indicate zones that move in the ranking.\n",
      "Zones ordered by mean rank. White columns are months with no acquisition.")
  ) +
  theme_pub(base_size = 9, grid = "none") +
  theme(legend.position = "right")

save_fig(p1, "fig01_rank_heatmap", DIR_OUT, 13, 12.5, TAG)


# --- Figure 2: mean rank with dispersion -------------------------------------
# A dot-and-whisker chart. The horizontal extent of each whisker is one
# standard deviation of the rank, so a short whisker is a persistent zone.

p2 <- persistence |>
  dplyr::arrange(mean_rank) |>
  dplyr::mutate(zone_name = factor(zone_name, levels = zone_name)) |>
  ggplot(aes(x = mean_rank, y = zone_name)) +
  geom_errorbarh(aes(xmin = mean_rank - sd_rank, xmax = mean_rank + sd_rank),
                 colour = "grey65", height = 0, linewidth = 0.5) +
  geom_point(aes(colour = lst_mean_all), size = 2.3) +
  scale_colour_gradientn(colours = PAL_LST, name = "Long-run\nmean (\u00b0C)") +
  labs(
    title    = sprintf("Mean position in the temperature ranking, %s", AOI_LABEL),
    subtitle = "Point is the mean rank, whisker is one standard deviation. Short whiskers identify zones whose position is stable across the record.",
    x = sprintf("Mean rank  (1 = coolest, %d = warmest)", n_zones), y = NULL,
    caption = "A zone can have a stable rank and still be unremarkable in absolute terms; the colour shows where it sits on the temperature scale."
  ) +
  theme_pub(base_size = 9, grid = "x") +
  theme(legend.position = "right")

save_fig(p2, "fig02_mean_rank_with_sd", DIR_OUT, 8.5, 12.5, TAG)


# --- Figure 3: time in the extreme groups ------------------------------------

pct_long <- persistence |>
  dplyr::select(zone_name, pct_top_warm, pct_top_cool) |>
  tidyr::pivot_longer(c(pct_top_warm, pct_top_cool),
                      names_to = "group", values_to = "pct") |>
  dplyr::mutate(
    group = dplyr::recode(group,
                          pct_top_warm = sprintf("Among the %d warmest", N_TOP),
                          pct_top_cool = sprintf("Among the %d coolest", N_TOP)),
    zone_name = stats::reorder(zone_name, pct, FUN = max))

p3 <- ggplot(pct_long, aes(x = pct, y = zone_name, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_vline(xintercept = 50, colour = "grey45", linetype = "dashed",
             linewidth = 0.4) +
  geom_vline(xintercept = 80, colour = "grey20", linetype = "dotted",
             linewidth = 0.5) +
  scale_fill_manual(
    values = stats::setNames(c("#c0392b", "#2c7bb6"),
                             c(sprintf("Among the %d warmest", N_TOP),
                               sprintf("Among the %d coolest", N_TOP))),
    name = NULL) +
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(
    title    = sprintf("Time spent in the temperature extremes, %s", AOI_LABEL),
    subtitle = sprintf("Share of observed months in which each zone was among the %d warmest or the %d coolest.",
                       N_TOP, N_TOP),
    x = "Share of observed months", y = NULL,
    caption = "Dashed line at 50%, dotted at 80%. Zones passing the dotted line are stable thermal extremes rather than occasional ones."
  ) +
  theme_pub(base_size = 9, grid = "x") +
  theme(legend.position = "bottom")

save_fig(p3, "fig03_time_in_extremes", DIR_OUT, 9, 12.5, TAG)


# --- Figure 4: mean anomaly relative to the city -----------------------------
# This is the covariate itself: the quantity actually fed to the model.

p4 <- persistence |>
  dplyr::arrange(delta_vs_city) |>
  dplyr::mutate(zone_name = factor(zone_name, levels = zone_name)) |>
  ggplot(aes(x = delta_vs_city, y = zone_name, fill = delta_vs_city)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 0, colour = "grey25", linewidth = 0.4) +
  scale_fill_gradient2(low = "#2c7bb6", mid = "grey95", high = "#c0392b",
                       midpoint = 0, name = "\u0394 (\u00b0C)") +
  labs(
    title    = sprintf("Mean surface temperature anomaly by zone, %s", AOI_LABEL),
    subtitle = "Difference between each zone and the city-wide mean of the same month, averaged over the record. This is the covariate used by the spatial model.",
    x = "Zone minus city mean (\u00b0C)", y = NULL,
    caption = "Positive values identify persistent heat islands. Because the city mean of the same acquisition is subtracted, any bias common to a whole scene cancels here."
  ) +
  theme_pub(base_size = 9, grid = "x") +
  theme(legend.position = "right")

save_fig(p4, "fig04_anomaly_vs_city", DIR_OUT, 8.5, 12.5, TAG)


# --- Figure 5: position against stability ------------------------------------
# A two-dimensional summary of the whole analysis. The horizontal axis is where
# a zone sits, the vertical axis is how much it moves. The interesting corners
# are the bottom two: stable extremes.

label_layer <- if (requireNamespace("ggrepel", quietly = TRUE)) {
  ggrepel::geom_text_repel(aes(label = zone_name), size = 2.5,
                           max.overlaps = 18, colour = "grey25",
                           segment.colour = "grey70", segment.size = 0.25,
                           min.segment.length = 0.2, box.padding = 0.25)
} else {
  geom_text(aes(label = zone_name), size = 2.5, colour = "grey25",
            check_overlap = TRUE, vjust = -0.8)
}

p5 <- ggplot(persistence, aes(x = mean_rank, y = sd_rank)) +
  geom_point(aes(colour = lst_mean_all), size = 3, alpha = 0.9) +
  label_layer +
  scale_colour_gradientn(colours = PAL_LST, name = "Long-run\nmean (\u00b0C)") +
  labs(
    title    = sprintf("Stability map of the temperature ranking, %s", AOI_LABEL),
    subtitle = "Bottom left: consistently cool. Bottom right: consistently warm. Upper region: zones whose position shifts from month to month.",
    x = sprintf("Mean rank  (1 = coolest, %d = warmest)", n_zones),
    y = "Standard deviation of the rank",
    caption = "Zones in the two lower corners are the ones a spatial exposure term describes well; those higher up carry more month-to-month noise."
  ) +
  theme_pub() +
  theme(legend.position = "right")

save_fig(p5, "fig05_stability_map", DIR_OUT, 10, 8, TAG)


# --- Figure 6: spatial spread over time --------------------------------------

p6 <- ggplot(amplitude, aes(x = date, y = amplitude)) +
  geom_line(colour = "grey70", linewidth = 0.4) +
  geom_point(aes(colour = amplitude), size = 1.9) +
  geom_hline(yintercept = mean(amplitude$amplitude), linetype = "dashed",
             colour = "grey30", linewidth = 0.4) +
  annotate("text", x = min(amplitude$date), y = mean(amplitude$amplitude),
           label = sprintf(" mean %.1f \u00b0C", mean(amplitude$amplitude)),
           hjust = 0, vjust = -0.6, size = 3, colour = "grey30") +
  scale_colour_gradientn(colours = PAL_LST, name = "Spread\n(\u00b0C)") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y",
               expand = expansion(mult = 0.01)) +
  labs(
    title    = sprintf("Spatial spread of surface temperature over time, %s", AOI_LABEL),
    subtitle = "Difference between the warmest and coolest zone on each acquisition date.",
    x = NULL, y = "Warmest minus coolest zone (\u00b0C)",
    caption = paste0(
      "This is the contrast a spatial model has available on a given date. ",
      "The spread widens in summer, when strong insolation amplifies the\n",
      "difference between vegetated and built-up surfaces.")
  ) +
  theme_pub(grid = "y") +
  theme(legend.position = "right")

save_fig(p6, "fig06_spatial_spread_over_time", DIR_OUT, 11, 5.5, TAG)


message(sprintf("%s Done: %s", TAG, Sys.time()))
message(sprintf("%s Outputs in %s\n", TAG, DIR_OUT))
