# =============================================================================
# MODULE 08 - CARTOGRAPHIC VISUALISATION OF NDVI AND LST
# =============================================================================
# Produces maps and charts at two spatial scales (pixel, zone) and two temporal
# scales (monthly, annual). Every label, title, legend and caption is in
# ENGLISH, so figures are publication-ready without editing.
#
# Outputs under OUTPUT_DIR/plots/:
#   pixel/monthly/  ndvi_pixel_YYYY_MM.png | lst_pixel_YYYY_MM.png
#   pixel/annual/   ndvi_pixel_YYYY.png    | lst_pixel_YYYY.png
#   zones/monthly/  ndvi_zone_YYYY_MM.png  | lst_zone_YYYY_MM.png
#   zones/annual/   ndvi_zone_YYYY.png     | lst_zone_YYYY.png
#                   ndvi_zone_annual_series.png | lst_zone_annual_series.png
#                   ndvi_zone_heatmap.png       | lst_zone_heatmap.png
#   timeseries/     aoi_monthly_series.png
#
# DEPENDENCIES - DELIBERATELY MINIMAL
#   This module needs only ggplot2, scales and (optionally) ggrepel.
#   Earlier versions required tidyterra and ggspatial; both were dropped
#   because they are frequently unavailable on cluster R installations and
#   would silently disable ALL figures. The north arrow and scale bar are
#   drawn directly with ggplot2 primitives instead, and rasters are converted
#   to data frames and drawn with geom_raster().
#
# ### THINGS YOU MAY WANT TO EDIT ARE MARKED ### EDIT ### BELOW. ###
# The most likely one is the COLOUR SCALE LIMITS (Section 1): they are fixed
# so that maps from different dates are directly comparable.
#
# Set GENERATE_PLOTS = FALSE in 00_config.R to skip this module entirely.
# =============================================================================

source("00_config.R")
source("04_clip_aoi.R")

# =============================================================================
# 0. PLOTTING PACKAGES (loaded lazily, so a missing graphics stack never
#    prevents the tables from being produced)
# =============================================================================

.HAS_GGREPEL <- FALSE

.ensure_plot_packages <- function() {

  required <- c("ggplot2", "scales")
  missing_pkgs <- required[!sapply(required, requireNamespace, quietly = TRUE)]

  if (length(missing_pkgs) > 0) {
    target_lib <- if (is.null(R_LIBS_USER)) .libPaths()[1] else R_LIBS_USER
    log_msg(sprintf("[08_plots] Installing: %s", paste(missing_pkgs, collapse = ", ")))
    try(install.packages(missing_pkgs, lib = target_lib,
                         repos = "https://cloud.r-project.org", quiet = TRUE),
        silent = TRUE)
  }

  ok <- all(sapply(required, requireNamespace, quietly = TRUE))
  if (!ok) {
    log_msg("[08_plots] ggplot2/scales unavailable - skipping all figures.", "WARN")
    return(FALSE)
  }
  suppressPackageStartupMessages({
    library(ggplot2, quietly = TRUE); library(scales, quietly = TRUE)
  })

  # ggrepel is optional: only used for zone labels
  .HAS_GGREPEL <<- requireNamespace("ggrepel", quietly = TRUE)
  TRUE
}

# =============================================================================
# 1. PALETTES, SCALE LIMITS AND THEME
# =============================================================================

# ### EDIT ### - colour ramps.
# NDVI: red (bare soil / water) -> yellow -> green (dense vegetation)
PAL_NDVI <- c("#d73027", "#fc8d59", "#fee090", "#d9ef8b", "#91cf60", "#1a9641")
# LST: blue (cool) -> white -> orange -> red (heat islands)
PAL_LST  <- c("#2b83ba", "#abdda4", "#ffffbf", "#fdae61", "#d7191c")

# ### EDIT ### - FIXED COLOUR-SCALE LIMITS.
#
# Deliberately FIXED rather than computed per map, so that any two maps in the
# series are visually comparable: the same colour always means the same value.
# The cost is that values outside the range are squished to the end colours.
#
# Choose them from your own data: run the pipeline once, open
# tables/monthly/ndvi_lst_monthly.csv, and pick limits covering roughly the
# 5th-95th percentile of the series.
#
#   Humid subtropical city (default) .. NDVI c(-0.1, 0.8)  LST c(15, 55)
#   Semi-arid / desert ................ NDVI c(-0.1, 0.5)  LST c(20, 70)
#   Temperate with winter ............. NDVI c(-0.1, 0.9)  LST c(-5, 45)
#   Dense forest ...................... NDVI c( 0.2, 0.95) LST c(15, 40)
# Set either to NULL for AUTOMATIC limits, computed from the 2nd-98th
# percentile of your own data on the first run. Automatic limits make a first
# run look right immediately, but they change if you later add scenes, so fix
# them to explicit numbers before producing final figures for publication.
PLOT_LIMITS_NDVI <- NULL          # e.g. c(-0.1, 0.8)
PLOT_LIMITS_LST  <- NULL          # e.g. c(15, 55)

# Filled in at run time when the above are NULL (do not edit).
.AUTO_LIMITS <- new.env(parent = emptyenv())

# ### EDIT ### - zone name labels on choropleth maps.
# Useful for a handful of zones, unreadable for hundreds.
PLOT_ZONE_LABELS    <- TRUE
PLOT_ZONE_LABEL_MAX <- 60     # skip labels above this many zones

# ### EDIT ### - output resolution. 300 dpi is the usual journal minimum.
PLOT_DPI    <- 300
PLOT_WIDTH  <- 9      # inches
PLOT_HEIGHT <- 8

#' Clean cartographic theme (map figures).
theme_carto <- function(base_size = 11) {
  ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = base_size + 2,
                                            hjust = 0.5,
                                            margin = ggplot2::margin(b = 4, t = 6)),
      plot.subtitle = ggplot2::element_text(size = base_size - 1, hjust = 0.5,
                                            colour = "grey40",
                                            margin = ggplot2::margin(b = 6)),
      plot.caption  = ggplot2::element_text(size = base_size - 3, colour = "grey55",
                                            hjust = 1,
                                            margin = ggplot2::margin(t = 4, b = 2)),
      plot.margin       = ggplot2::margin(8, 8, 6, 8),
      legend.position   = "right",
      legend.title      = ggplot2::element_text(face = "bold", size = base_size - 1),
      legend.text       = ggplot2::element_text(size = base_size - 2),
      legend.key.width  = ggplot2::unit(0.4, "cm"),
      legend.key.height = ggplot2::unit(1.6, "cm"),
      panel.border      = ggplot2::element_rect(colour = "grey60", fill = NA,
                                                linewidth = 0.4),
      plot.background   = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

#' Clean theme for non-map charts (time series, heatmaps).
theme_chart <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold"),
      plot.subtitle    = ggplot2::element_text(colour = "grey40"),
      plot.caption     = ggplot2::element_text(size = base_size - 4, colour = "grey55",
                                               hjust = 1),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      legend.position  = "right"
    )
}

# =============================================================================
# 2. INTERNAL HELPERS
# =============================================================================

.create_plot_dirs <- function() {
  dirs <- c(
    file.path(OUTPUT_DIR, "plots", "pixel", "monthly"),
    file.path(OUTPUT_DIR, "plots", "pixel", "annual"),
    file.path(OUTPUT_DIR, "plots", "zones", "monthly"),
    file.path(OUTPUT_DIR, "plots", "zones", "annual"),
    file.path(OUTPUT_DIR, "plots", "timeseries")
  )
  lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  invisible(dirs)
}

#' Build a title/filename pair for a period. Month names are English because
#' the pipeline sets no locale-specific formatting for figures.
label_period <- function(year, month = NA, month_name = NA) {
  if (!is.na(month)) {
    nm <- month.name[month]     # always English
    list(title = sprintf("%s %d", nm, year), file = sprintf("%04d_%02d", year, month))
  } else {
    list(title = as.character(year), file = sprintf("%04d", year))
  }
}

save_plot <- function(p, path, width = PLOT_WIDTH, height = PLOT_HEIGHT,
                      dpi = PLOT_DPI) {
  ok <- tryCatch({
    ggplot2::ggsave(path, plot = p, width = width, height = height,
                    dpi = dpi, units = "in", bg = "white")
    TRUE
  }, error = function(e) {
    log_msg(sprintf("[08_plots] ggsave failed for %s: %s",
                    basename(path), conditionMessage(e)), "ERROR")
    FALSE
  })
  if (ok) log_msg(sprintf("[08_plots] Saved: %s", basename(path)))
  invisible(path)
}

#' Caption line shared by every figure (provenance + date).
.map_caption <- function() {
  sprintf("Landsat Collection 2 Level-2 (USGS)  |  CRS: %s  |  Generated %s",
          CRS_TARGET, format(Sys.Date(), "%Y-%m-%d"))
}

#' Variable-specific plotting parameters (English labels throughout).
.var_params <- function(variable) {
  if (toupper(variable) == "NDVI") {
    lim <- if (!is.null(PLOT_LIMITS_NDVI)) PLOT_LIMITS_NDVI
           else get0("ndvi", envir = .AUTO_LIMITS, ifnotfound = c(-0.1, 0.8))
    brk <- pretty(lim, 6)
    list(pal = PAL_NDVI, lim = lim, leg_title = "NDVI",
         breaks = brk, labels = sprintf("%.2f", brk),
         axis_title = "NDVI (dimensionless)", unit = "")
  } else {
    lim <- if (!is.null(PLOT_LIMITS_LST)) PLOT_LIMITS_LST
           else get0("lst", envir = .AUTO_LIMITS, ifnotfound = c(15, 55))
    brk <- pretty(lim, 5)
    list(pal = PAL_LST, lim = lim, leg_title = "LST (\u00b0C)",
         breaks = brk, labels = paste0(brk, "\u00b0"),
         axis_title = "Land surface temperature (\u00b0C)", unit = " \u00b0C")
  }
}

#' Per-period pixelwise mean of the saved rasters (display only).
mean_rasters_period <- function(df_scenes, type, dir_ras) {
  suffix  <- if (type == "ndvi") "_NDVI.tif" else "_LST_Celsius.tif"
  rasters <- list()
  for (sid in df_scenes$scene_id) {
    f <- file.path(dir_ras, paste0(sid, suffix))
    if (file.exists(f)) rasters <- c(rasters, list(terra::rast(f)))
  }
  if (length(rasters) == 0) return(NULL)

  template <- rasters[[1]]
  for (k in seq_along(rasters))
    if (!terra::compareGeom(rasters[[k]], template, stopOnError = FALSE))
      rasters[[k]] <- terra::resample(rasters[[k]], template, method = "bilinear")

  m <- terra::app(terra::rast(rasters), mean, na.rm = TRUE)
  names(m) <- type
  m
}

#' Convert a SpatRaster to a data frame ggplot2 can draw with geom_raster().
#' (Replaces tidyterra::geom_spatraster, removing that dependency.)
.raster_to_df <- function(r) {
  df <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  if (ncol(df) < 3 || nrow(df) == 0) return(NULL)
  names(df)[3] <- "value"
  df
}

#' Convert an sf polygon layer to a data frame of ring coordinates, so that
#' boundaries can be drawn with geom_polygon() without needing geom_sf.
.sf_to_path_df <- function(x) {
  geo <- sf::st_geometry(x)
  out <- list(); k <- 0
  for (i in seq_along(geo)) {
    g <- sf::st_cast(sf::st_sfc(geo[[i]], crs = sf::st_crs(x)), "POLYGON",
                     warn = FALSE)
    for (j in seq_along(g)) {
      rings <- sf::st_coordinates(g[j])
      for (rg in unique(rings[, "L1"])) {
        k <- k + 1
        sub <- rings[rings[, "L1"] == rg, , drop = FALSE]
        out[[k]] <- data.frame(x = sub[, "X"], y = sub[, "Y"],
                               grp = sprintf("%d_%d_%d", i, j, rg),
                               feature = i)
      }
    }
  }
  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}

#' Draw a north arrow and a graphic scale bar directly with ggplot2 layers.
#' (Replaces ggspatial::annotation_north_arrow / annotation_scale.)
#'
#' @param xr,yr numeric range of the plotting extent, in CRS units (metres)
.add_map_furniture <- function(p, xr, yr) {

  dx <- diff(xr); dy <- diff(yr)

  # --- north arrow: triangle + "N", top-right ---
  ax <- xr[2] - 0.05 * dx
  ay <- yr[2] - 0.06 * dy
  as <- 0.030 * dy                      # arrow half-height
  aw <- 0.011 * dx                      # arrow half-width

  arrow_df <- data.frame(
    x = c(ax, ax - aw, ax, ax + aw),
    y = c(ay + as, ay - as, ay - as * 0.45, ay - as)
  )

  # --- scale bar: "nice" round length, bottom-left ---
  target <- dx * 0.25
  nice   <- c(10, 20, 25, 50, 100, 200, 250, 500,
              1000, 2000, 2500, 5000, 10000, 20000, 25000, 50000,
              100000, 200000, 500000)
  bar_m  <- nice[which.min(abs(nice - target))]
  bx0    <- xr[1] + 0.05 * dx
  by     <- yr[1] + 0.05 * dy
  bh     <- 0.010 * dy

  bar_lab <- if (bar_m >= 1000) sprintf("%g km", bar_m / 1000) else sprintf("%g m", bar_m)

  # two-tone bar (dark | white) for legibility over any background
  seg <- data.frame(
    xmin = c(bx0, bx0 + bar_m / 2),
    xmax = c(bx0 + bar_m / 2, bx0 + bar_m),
    ymin = by, ymax = by + bh,
    fill_col = c("grey20", "white")
  )

  p +
    ggplot2::geom_polygon(data = arrow_df,
                          ggplot2::aes(x = x, y = y),
                          fill = "grey20", colour = "grey20",
                          linewidth = 0.3, inherit.aes = FALSE) +
    ggplot2::annotate("text", x = ax, y = ay + as * 1.55, label = "N",
                      size = 3.4, fontface = "bold", colour = "grey20") +
    ggplot2::geom_rect(data = seg,
                       ggplot2::aes(xmin = xmin, xmax = xmax,
                                    ymin = ymin, ymax = ymax),
                       fill = seg$fill_col, colour = "grey20",
                       linewidth = 0.3, inherit.aes = FALSE) +
    ggplot2::annotate("text", x = bx0 + bar_m / 2, y = by + bh * 2.6,
                      label = bar_lab, size = 2.9, colour = "grey20")
}

# =============================================================================
# 3. PIXEL-SCALE MAP (30 m)
# =============================================================================

plot_raster_period <- function(r, variable, period, boundary_sf, dir_out) {

  if (is.null(r)) {
    log_msg(sprintf("[08_plots] NULL raster for %s %s - skipping.",
                    variable, period$file), "WARN")
    return(invisible(NULL))
  }

  df <- .raster_to_df(r)
  if (is.null(df)) {
    log_msg(sprintf("[08_plots] No valid pixels for %s %s - skipping.",
                    variable, period$file), "WARN")
    return(invisible(NULL))
  }

  vp   <- .var_params(variable)
  bdf  <- .sf_to_path_df(boundary_sf)
  xr   <- range(df$x); yr <- range(df$y)

  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = df, ggplot2::aes(x = x, y = y, fill = value))

  if (!is.null(bdf))
    p <- p + ggplot2::geom_polygon(data = bdf,
                                   ggplot2::aes(x = x, y = y, group = grp),
                                   fill = NA, colour = "grey15", linewidth = 0.5)

  p <- p +
    ggplot2::scale_fill_gradientn(
      colours  = vp$pal, limits = vp$lim,
      oob      = scales::squish,     # out-of-range -> end colour, not blank
      name     = vp$leg_title, breaks = vp$breaks, labels = vp$labels,
      na.value = "transparent",
      guide    = ggplot2::guide_colourbar(ticks.colour = "grey30",
                                          frame.colour = "grey30")) +
    ggplot2::coord_equal(expand = FALSE)

  p <- .add_map_furniture(p, xr, yr)

  p <- p +
    ggplot2::labs(
      title    = sprintf("%s \u2014 %s, %s", vp$leg_title, AOI_NAME, period$title),
      subtitle = sprintf("Pixel-level mean (30 \u00d7 30 m)  |  fixed colour scale [%.1f, %.1f]%s",
                         vp$lim[1], vp$lim[2], vp$unit),
      caption  = .map_caption()) +
    theme_carto()

  save_plot(p, file.path(dir_out, sprintf("%s_pixel_%s.png",
                                          tolower(variable), period$file)))
}

# =============================================================================
# 4. CHOROPLETH MAP (ZONE SCALE)
# =============================================================================

plot_zone_period <- function(df_zone, zones_sf, variable, col_val,
                             period, dir_out, labels = FALSE) {

  vp <- .var_params(variable)

  vals <- df_zone[, c("zone_name", col_val)]
  names(vals)[2] <- "value"
  map_sf <- dplyr::left_join(zones_sf, vals, by = "zone_name")

  if (all(is.na(map_sf$value))) {
    log_msg(sprintf("[08_plots] No zone data for %s %s - skipping.",
                    variable, period$file), "WARN")
    return(invisible(NULL))
  }

  pdf_ <- .sf_to_path_df(map_sf)
  if (is.null(pdf_)) return(invisible(NULL))
  pdf_$value <- map_sf$value[pdf_$feature]

  xr <- range(pdf_$x); yr <- range(pdf_$y)

  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = pdf_,
                          ggplot2::aes(x = x, y = y, group = grp, fill = value),
                          colour = "grey30", linewidth = 0.25) +
    ggplot2::scale_fill_gradientn(
      colours  = vp$pal, limits = vp$lim, oob = scales::squish,
      name     = vp$leg_title, breaks = vp$breaks, labels = vp$labels,
      na.value = "grey88",            # grey = zone with no valid data
      guide    = ggplot2::guide_colourbar(ticks.colour = "grey30",
                                          frame.colour = "grey30")) +
    ggplot2::coord_equal(expand = FALSE)

  if (isTRUE(labels) && isTRUE(PLOT_ZONE_LABELS) &&
      nrow(map_sf) <= PLOT_ZONE_LABEL_MAX) {
    cent <- suppressWarnings(sf::st_coordinates(sf::st_point_on_surface(
      sf::st_geometry(map_sf))))
    lab_df <- data.frame(x = cent[, 1], y = cent[, 2],
                         label = map_sf$zone_name)
    if (.HAS_GGREPEL) {
      p <- p + ggrepel::geom_text_repel(
        data = lab_df, ggplot2::aes(x = x, y = y, label = label),
        size = 2.4, colour = "grey10", segment.size = 0.2,
        min.segment.length = 0.2, max.overlaps = 30, inherit.aes = FALSE)
    } else {
      p <- p + ggplot2::geom_text(
        data = lab_df, ggplot2::aes(x = x, y = y, label = label),
        size = 2.4, colour = "grey10", inherit.aes = FALSE)
    }
  }

  p <- .add_map_furniture(p, xr, yr)

  p <- p +
    ggplot2::labs(
      title    = sprintf("%s by %s \u2014 %s, %s",
                         vp$leg_title, ZONE_LABEL, AOI_NAME, period$title),
      subtitle = sprintf("Pixel-count weighted mean per %s  |  grey = no valid data",
                         ZONE_LABEL),
      caption  = .map_caption()) +
    theme_carto()

  save_plot(p, file.path(dir_out, sprintf("%s_zone_%s.png",
                                          tolower(variable), period$file)))
}

# =============================================================================
# 5. ANNUAL TIME SERIES BY ZONE
# =============================================================================

plot_zone_series <- function(df_annual, variable, col_val, dir_out) {

  vp <- .var_params(variable)

  df <- df_annual[!is.na(df_annual[[col_val]]), ]
  if (nrow(df) == 0) return(invisible(NULL))
  df$value <- df[[col_val]]

  aoi_mean <- df |>
    dplyr::group_by(ano) |>
    dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = ano, y = value, group = zone_name)) +
    ggplot2::geom_line(colour = "grey70", linewidth = 0.4, alpha = 0.9) +
    ggplot2::geom_point(colour = "grey60", size = 0.8) +
    ggplot2::geom_line(data = aoi_mean,
                       ggplot2::aes(x = ano, y = value, group = 1),
                       colour = "#b2182b", linewidth = 1.2, inherit.aes = FALSE) +
    ggplot2::geom_point(data = aoi_mean,
                        ggplot2::aes(x = ano, y = value, group = 1),
                        colour = "#b2182b", size = 2, inherit.aes = FALSE) +
    ggplot2::scale_x_continuous(breaks = sort(unique(df$ano))) +
    ggplot2::labs(
      title    = sprintf("Annual %s by %s \u2014 %s", vp$leg_title, ZONE_LABEL, AOI_NAME),
      subtitle = sprintf("Grey lines: individual %ss.  Red line: %s-wide mean.",
                         ZONE_LABEL, AOI_NAME),
      x = "Year", y = vp$axis_title,
      caption = .map_caption()) +
    theme_chart()

  save_plot(p, file.path(dir_out, sprintf("%s_zone_annual_series.png",
                                          tolower(variable))),
            width = 10, height = 5.5)
}

# =============================================================================
# 6. ZONE x TIME HEATMAP
# =============================================================================

#' Compact overview of the whole series: one row per zone, one column per
#' month, colour = value. Reveals data gaps and zone-level anomalies at a
#' glance, which the line chart hides.
plot_zone_heatmap <- function(df_monthly, variable, col_val, dir_out) {

  vp <- .var_params(variable)
  df <- df_monthly
  df$value <- df[[col_val]]
  df <- df[!is.na(df$zone_name), ]
  if (all(is.na(df$value))) return(invisible(NULL))

  df$date <- as.Date(sprintf("%d-%02d-15", df$ano, df$mes))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = date, y = zone_name, fill = value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.15) +
    ggplot2::scale_fill_gradientn(
      colours = vp$pal, limits = vp$lim, oob = scales::squish,
      name = vp$leg_title, breaks = vp$breaks, labels = vp$labels,
      na.value = "grey92") +
    ggplot2::scale_x_date(date_labels = "%Y-%m", date_breaks = "3 months") +
    ggplot2::labs(
      title    = sprintf("Monthly %s by %s \u2014 %s", vp$leg_title, ZONE_LABEL, AOI_NAME),
      subtitle = "Grey cells: no cloud-free observation available in that month",
      x = "Month", y = NULL, caption = .map_caption()) +
    theme_chart() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                                       size = 7),
                   axis.text.y = ggplot2::element_text(size = 7),
                   panel.grid  = ggplot2::element_blank())

  save_plot(p, file.path(dir_out, sprintf("%s_zone_heatmap.png", tolower(variable))),
            width = 11, height = max(3.5, 0.22 * length(unique(df$zone_name)) + 2.5))
}

# =============================================================================
# 7. AOI-WIDE MONTHLY SERIES (NDVI + LST, dual panel)
# =============================================================================

plot_aoi_series <- function(df_monthly, dir_out) {

  df <- df_monthly
  df$date <- as.Date(sprintf("%d-%02d-15", df$ano, df$mes))
  if (nrow(df) == 0 || all(is.na(df$ndvi_media) & is.na(df$lst_c_media)))
    return(invisible(NULL))

  # NOTE: NA rows are deliberately KEPT. geom_line() breaks the line at NA,
  # so months without a cloud-free observation appear as genuine gaps rather
  # than being silently interpolated across - which is what the subtitle
  # promises the reader.
  long <- rbind(
    data.frame(date = df$date, value = df$ndvi_media,
               panel = "NDVI (dimensionless)"),
    data.frame(date = df$date, value = df$lst_c_media,
               panel = "Land surface temperature (\u00b0C)")
  )

  # Two layers, deliberately distinguished:
  #   - a faint DASHED connector through all observations, so the eye can
  #     follow the series even when months are missing;
  #   - a solid SOLID line that breaks at NA, so real gaps stay visible;
  #   - points marking the months that actually have an observation.
  # The reader can therefore never mistake an interpolated segment for data.
  obs <- long[!is.na(long$value), ]

  p <- ggplot2::ggplot(long, ggplot2::aes(x = date, y = value)) +
    ggplot2::geom_line(data = obs, colour = "grey70", linewidth = 0.4,
                       linetype = "22", na.rm = TRUE) +
    ggplot2::geom_line(colour = "#2166ac", linewidth = 0.7, na.rm = TRUE) +
    ggplot2::geom_point(data = obs, colour = "#2166ac", size = 1.5, na.rm = TRUE) +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_y") +
    ggplot2::scale_x_date(date_labels = "%Y-%m", date_breaks = "3 months") +
    ggplot2::labs(
      title    = sprintf("Monthly NDVI and land surface temperature \u2014 %s", AOI_NAME),
      subtitle = "Area-wide weighted mean  \u2014  points: observed months  \u2014  dashed: interpolated across gaps",
      x = "Month", y = NULL, caption = .map_caption()) +
    theme_chart() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                                       size = 8),
                   strip.text  = ggplot2::element_text(face = "bold"))

  save_plot(p, file.path(dir_out, "aoi_monthly_series.png"),
            width = 10, height = 6.5)
}

# =============================================================================
# 8. BATCH DRIVERS
# =============================================================================

plot_all_pixel <- function(dir_ndvi, dir_lst, df_scenes, boundary_sf,
                           scale = "monthly") {

  dir_out <- file.path(OUTPUT_DIR, "plots", "pixel", scale)

  groups <- if (scale == "monthly")
    split(df_scenes, list(df_scenes$ano, df_scenes$mes), drop = TRUE)
  else
    split(df_scenes, df_scenes$ano, drop = TRUE)

  for (g in groups) {
    if (nrow(g) == 0) next
    period <- if (scale == "monthly")
      label_period(g$ano[1], g$mes[1]) else label_period(g$ano[1])

    plot_raster_period(mean_rasters_period(g, "ndvi", dir_ndvi),
                       "NDVI", period, boundary_sf, dir_out)
    plot_raster_period(mean_rasters_period(g, "lst", dir_lst),
                       "LST", period, boundary_sf, dir_out)
    gc(verbose = FALSE)
  }
}

plot_all_zone <- function(df_zone, zones_sf, scale = "monthly") {

  dir_out <- file.path(OUTPUT_DIR, "plots", "zones", scale)

  if (scale == "monthly") {
    df <- df_zone[!is.na(df_zone$ndvi_media) | !is.na(df_zone$lst_c_media), ]
    if (nrow(df) == 0) return(invisible(NULL))
    for (g in split(df, list(df$ano, df$mes), drop = TRUE)) {
      period <- label_period(g$ano[1], g$mes[1])
      plot_zone_period(g, zones_sf, "NDVI", "ndvi_media",  period, dir_out, labels = FALSE)
      plot_zone_period(g, zones_sf, "LST",  "lst_c_media", period, dir_out, labels = FALSE)
    }
    # heatmaps summarise the whole monthly series in two figures
    plot_zone_heatmap(df_zone, "NDVI", "ndvi_media",  dir_out)
    plot_zone_heatmap(df_zone, "LST",  "lst_c_media", dir_out)
  } else {
    df <- df_zone[!is.na(df_zone$ndvi_media_anual) | !is.na(df_zone$lst_media_anual), ]
    if (nrow(df) == 0) return(invisible(NULL))
    for (g in split(df, df$ano, drop = TRUE)) {
      period <- label_period(g$ano[1])
      plot_zone_period(g, zones_sf, "NDVI", "ndvi_media_anual", period, dir_out, labels = TRUE)
      plot_zone_period(g, zones_sf, "LST",  "lst_media_anual",  period, dir_out, labels = TRUE)
    }
    plot_zone_series(df, "NDVI", "ndvi_media_anual", dir_out)
    plot_zone_series(df, "LST",  "lst_media_anual",  dir_out)
  }
}

# =============================================================================
# ORCHESTRATOR
# =============================================================================

#' Generate every map and chart.
#'
#' @param df_scene_stats  per-scene statistics (module 05)
#' @param df_zone_monthly monthly per-zone table (module 07)
#' @param df_zone_annual  annual per-zone table (module 07)
#' @param zones_sf        sf zones from load_aoi_zones()
#' @param df_monthly      AOI-wide monthly table (module 05), optional
run_plots <- function(df_scene_stats, df_zone_monthly, df_zone_annual, zones_sf,
                      df_monthly = NULL) {

  log_msg("[08_plots] === START - Figure generation ===")

  if (!.ensure_plot_packages()) {
    log_msg("[08_plots] Aborting figures; tables are unaffected.", "WARN")
    return(invisible(NULL))
  }

  .create_plot_dirs()

  # --- automatic colour limits, when PLOT_LIMITS_* are NULL -----------------
  # Uses the 2nd-98th percentile of the per-scene means, rounded outwards, so
  # that the scale covers the bulk of the series without being dragged by a
  # single extreme scene.
  .set_auto <- function(key, vals, digits) {
    vals <- vals[is.finite(vals)]
    if (length(vals) < 2) return(invisible(NULL))
    q <- stats::quantile(vals, c(0.02, 0.98), na.rm = TRUE)
    pad <- max(diff(q) * 0.15, 10^(-digits))
    lim <- c(floor((q[1] - pad) * 10^digits) / 10^digits,
             ceiling((q[2] + pad) * 10^digits) / 10^digits)
    assign(key, as.numeric(lim), envir = .AUTO_LIMITS)
    log_msg(sprintf("[08_plots] Automatic %s colour limits: [%.2f, %.2f]. Fix them in 08_plots.R for final figures.",
                    toupper(key), lim[1], lim[2]))
  }
  if (is.null(PLOT_LIMITS_NDVI) && !is.null(df_scene_stats))
    .set_auto("ndvi", df_scene_stats$ndvi_media, 2)
  if (is.null(PLOT_LIMITS_LST) && !is.null(df_scene_stats))
    .set_auto("lst", c(df_scene_stats$lst_c_media, df_scene_stats$lst_c_max), 0)

  boundary_sf <- load_aoi_boundary()
  dir_ndvi    <- file.path(OUTPUT_DIR, "rasters", "ndvi")
  dir_lst     <- file.path(OUTPUT_DIR, "rasters", "lst")

  # --- AOI-wide series ---
  if (!is.null(df_monthly)) {
    tryCatch(plot_aoi_series(df_monthly, file.path(OUTPUT_DIR, "plots", "timeseries")),
             error = function(e)
               log_msg(sprintf("[08_plots] ERROR in AOI series: %s", e$message), "ERROR"))
  }

  # --- pixel maps (need the saved GeoTIFFs) ---
  if (isTRUE(SAVE_CLIPPED_RASTERS)) {
    log_msg("[08_plots] --- Monthly pixel maps ---")
    tryCatch(plot_all_pixel(dir_ndvi, dir_lst, df_scene_stats, boundary_sf, "monthly"),
             error = function(e) log_msg(sprintf("[08_plots] ERROR monthly pixel: %s", e$message), "ERROR"))
    log_msg("[08_plots] --- Annual pixel maps ---")
    tryCatch(plot_all_pixel(dir_ndvi, dir_lst, df_scene_stats, boundary_sf, "annual"),
             error = function(e) log_msg(sprintf("[08_plots] ERROR annual pixel: %s", e$message), "ERROR"))
  } else {
    log_msg("[08_plots] Pixel maps skipped (SAVE_CLIPPED_RASTERS = FALSE).", "WARN")
  }

  # --- zone maps ---
  if (!is.null(df_zone_monthly)) {
    log_msg("[08_plots] --- Monthly choropleths + heatmaps ---")
    tryCatch(plot_all_zone(df_zone_monthly, zones_sf, "monthly"),
             error = function(e) log_msg(sprintf("[08_plots] ERROR monthly zone: %s", e$message), "ERROR"))
  }
  if (!is.null(df_zone_annual)) {
    log_msg("[08_plots] --- Annual choropleths + series ---")
    tryCatch(plot_all_zone(df_zone_annual, zones_sf, "annual"),
             error = function(e) log_msg(sprintf("[08_plots] ERROR annual zone: %s", e$message), "ERROR"))
  }

  log_msg("[08_plots] === END - Figure generation ===")
  invisible(NULL)
}
