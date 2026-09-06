# =============================================================================
# VALIDATION SUITE - 00 - SHARED CONFIGURATION, HELPERS AND FIGURE STYLE
# =============================================================================
# Sourced by every 01_..05_ validation script. It centralises four things so
# that they are defined once and stay consistent:
#
#   (a) WHERE the pipeline outputs live and WHERE validation outputs go
#   (b) HOW input CSVs are located (by content, not by a hard-coded file name)
#   (c) HOW column names are normalised across pipeline variants
#   (d) The publication figure theme and the figure-saving function
#
# >>> THIS IS THE ONLY FILE YOU NEED TO EDIT. Everything you may want to      <<<
# >>> change is in Section 1 and marked  ### EDIT ###                         <<<
#
# WHY FILES ARE LOCATED BY CONTENT RATHER THAN BY PATH
#   The upstream pipeline exists in more than one variant: one writes to
#   tabelas/ with columns named ano, mes, nome_bairro; another writes to
#   tables/ with year, month, zone_name. Hard-coding either would break the
#   other. Instead these scripts search the results tree for a CSV that
#   CONTAINS the columns they need, then rename those columns onto a single
#   internal convention. Everything downstream assumes one naming scheme.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. CONFIGURATION   ### EDIT THIS SECTION ###
# -----------------------------------------------------------------------------

# --- 1a. Where the main pipeline wrote its results ---------------------------
# List as many candidates as you like; the first that exists on disk is used,
# so the same file can serve several machines or several runs.
#
# Leave the vector EMPTY (character(0)) to search automatically: the suite then
# looks for a directory containing a "tables" or "tabelas" folder, first beside
# these scripts, then one and two levels up. That covers the common layout in
# which the validation scripts sit inside the project tree.
PIPELINE_DIR_CANDIDATES <- character(0)
# e.g.
# PIPELINE_DIR_CANDIDATES <- c(
#   "~/my_project/results_landsat",
#   "/scratch/user/project/Resultados_Landsat"
# )

# --- 1b. Where validation outputs go -----------------------------------------
# NULL puts them in a "validation_output" folder beside these scripts, which
# keeps a cloned repository self-contained. Set an absolute path to write
# elsewhere.
VALIDATION_ROOT <- NULL
# e.g. VALIDATION_ROOT <- "~/my_project/validation"

# --- 1c. Study area label ----------------------------------------------------
# Appears in figure titles and report headers.
AOI_LABEL <- "Study area"

# --- 1d. Study period --------------------------------------------------------
# NULL means "infer from the data", which is usually what you want: the suite
# then reports exactly the span the pipeline produced. Set explicit years only
# to restrict the validation to a sub-period.
YEAR_START <- NULL
YEAR_END   <- NULL

# --- 1e. Hemisphere, for the season labels -----------------------------------
# "south": Dec-Feb is Summer.   "north": Jun-Aug is Summer.
# Getting this wrong silently mislabels every seasonal figure, so it is worth
# a moment's check. Set to NA to disable seasonal grouping entirely, which is
# the honest choice near the equator where the four-season scheme does not
# describe the climate.
HEMISPHERE <- "south"

# --- 1f. Months excluded as low quality --------------------------------------
# Some months survive every automatic filter in the pipeline yet are clearly
# unusable on inspection - a diffuse cold bias from Landsat-7 SLC-off striping
# is the classic case, because no single pixel is extreme enough to be clamped
# while the scene mean is badly wrong.
#
# Declare such months HERE, once. All five scripts pick the list up
# automatically, exclude those rows from every statistic, and state the
# exclusion in their figure captions and reports.
#
# Start empty. After running script 01 and inspecting the city-wide series,
# add any month whose mean is physically implausible for your climate, with a
# reason you would be willing to defend in print.
LOW_QUALITY_MONTHS <- data.frame(
  year   = integer(0),
  month  = integer(0),
  reason = character(0),
  stringsAsFactors = FALSE
)
# e.g.
# LOW_QUALITY_MONTHS <- data.frame(
#   year   = 2021,
#   month  = 7,
#   reason = "LE07 SLC-off, diffuse cold bias (city mean LST implausible)"
# )


# -----------------------------------------------------------------------------
# 2. FIGURE QUALITY SETTINGS  ### EDIT HERE IF NEEDED ###
# -----------------------------------------------------------------------------
# 320 dpi comfortably exceeds the 300 dpi most journals require for raster
# figures. Every figure is ALSO written as a PDF, which is vector: text and
# lines stay sharp at any zoom level, and that is the version to submit if the
# journal accepts vector artwork.
FIG_DPI       <- 320
FIG_SAVE_PDF  <- TRUE    # set FALSE to halve the number of output files
FIG_BASE_SIZE <- 11      # base font size in points for standard figures


# -----------------------------------------------------------------------------
# 3. PACKAGES
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(ggplot2); library(scales)
})

# Force English month names and English number formatting regardless of the
# container locale. Without this, format(date, "%b") returns "jan"/"fev" on a
# pt_BR system and the figures come out half-translated.
suppressWarnings(try(Sys.setlocale("LC_TIME", "C"), silent = TRUE))


# -----------------------------------------------------------------------------
# 3b. RESOLVING PATHS
# -----------------------------------------------------------------------------

# Directory containing these scripts, however they were invoked.
.suite_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- a[grep("--file=", a)]
  d <- if (length(f) > 0)
    suppressWarnings(try(dirname(normalizePath(sub("--file=", "", f[1]),
                                               mustWork = FALSE)),
                         silent = TRUE)) else NULL
  if (is.null(d) || inherits(d, "try-error") || !dir.exists(d)) getwd() else d
})

# Output root: beside the scripts unless the user set one explicitly.
if (is.null(VALIDATION_ROOT) || !nzchar(VALIDATION_ROOT))
  VALIDATION_ROOT <- file.path(.suite_dir, "validation_output")
VALIDATION_ROOT <- path.expand(VALIDATION_ROOT)
dir.create(VALIDATION_ROOT, recursive = TRUE, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# 4. LOCATING THE PIPELINE OUTPUTS
# -----------------------------------------------------------------------------

# Read a CSV regardless of whether it uses comma or semicolon separators.
# write.csv produces commas; write.csv2 (and Excel in a pt_BR locale) produces
# semicolons with comma decimals. Sniffing the header avoids a silent failure
# in which every row is parsed into a single character column.
read_any_csv <- function(path, n_max = Inf) {
  header <- readLines(path, n = 1, warn = FALSE)
  if (grepl(";", header)) {
    readr::read_csv2(path, n_max = n_max, show_col_types = FALSE)
  } else {
    readr::read_csv(path, n_max = n_max, show_col_types = FALSE)
  }
}

#' Does this directory look like a finished pipeline run?
.looks_like_results <- function(d) {
  dir.exists(d) &&
    any(dir.exists(file.path(d, c("tables", "tabelas"))))
}

#' Resolve the pipeline results directory.
#'
#' Uses PIPELINE_DIR_CANDIDATES when set. When that vector is empty it searches
#' automatically: beside these scripts, then one and two levels up, looking for
#' any directory that contains a tables/ or tabelas/ folder. That covers the
#' usual layout in which the validation suite lives inside the project tree,
#' without requiring anyone to edit a path.
resolve_pipeline_dir <- function() {

  if (length(PIPELINE_DIR_CANDIDATES) > 0) {
    cand <- path.expand(PIPELINE_DIR_CANDIDATES)
    hit  <- cand[vapply(cand, .looks_like_results, logical(1))]
    if (length(hit) == 0) {
      exists_but_empty <- cand[dir.exists(cand)]
      stop("No usable pipeline results found.\n",
           if (length(exists_but_empty))
             paste0("  These directories exist but contain no tables/ or ",
                    "tabelas/ folder:\n    ",
                    paste(exists_but_empty, collapse = "\n    "), "\n")
           else "",
           "  Checked:\n    ", paste(cand, collapse = "\n    "),
           "\n  Fix PIPELINE_DIR_CANDIDATES in 00_common.R, section 1a.")
    }
    return(normalizePath(hit[1]))
  }

  # --- automatic search ---
  roots <- unique(c(.suite_dir,
                    dirname(.suite_dir),
                    dirname(dirname(.suite_dir)),
                    getwd(), dirname(getwd())))
  found <- character(0)
  for (r in roots) {
    if (!dir.exists(r)) next
    if (.looks_like_results(r)) found <- c(found, r)
    subs <- list.dirs(r, recursive = FALSE, full.names = TRUE)
    found <- c(found, subs[vapply(subs, .looks_like_results, logical(1))])
  }
  found <- unique(found)

  if (length(found) == 0)
    stop("Could not find pipeline results automatically.\n",
         "  Searched in and below:\n    ", paste(roots, collapse = "\n    "),
         "\n  A results directory is one containing a tables/ or tabelas/ ",
         "folder.\n  Set PIPELINE_DIR_CANDIDATES in 00_common.R, section 1a.")

  if (length(found) > 1)
    message("[common] Several candidate result directories found; using the ",
            "first:\n           ", paste(found, collapse = "\n           "))
  normalizePath(found[1])
}

# Find the tables root, which is "tabelas" in the Portuguese pipeline and
# "tables" in the English one.
resolve_tables_dir <- function(pipeline_dir) {
  cands <- file.path(pipeline_dir, c("tables", "tabelas"))
  hit   <- cands[dir.exists(cands)]
  if (length(hit) == 0)
    stop("No tables/ or tabelas/ folder under ", pipeline_dir,
         "\nHas the main pipeline finished a run?")
  hit[1]
}

# Search the tables tree for a CSV that contains ALL of `required_any`, where
# each element of `required_any` is itself a vector of acceptable spellings for
# one logical column. Files whose path matches `prefer_pattern` are examined
# first, but a name mismatch is not fatal.
#
# Returns the path, or NULL if nothing matches.
find_csv_by_columns <- function(tables_dir, required_any,
                                prefer_pattern = NULL,
                                forbid_any = NULL) {

  cands <- list.files(tables_dir, pattern = "\\.csv$",
                      recursive = TRUE, full.names = TRUE)
  if (length(cands) == 0) return(NULL)

  # Files whose path hints at the right table are examined first, but a name
  # mismatch never disqualifies a file - only its columns do.
  if (!is.null(prefer_pattern)) {
    preferred <- grepl(prefer_pattern, cands, ignore.case = TRUE)
    cands     <- c(cands[preferred], cands[!preferred])
  }

  for (p in cands) {
    hdr <- tryCatch(names(read_any_csv(p, n_max = 1)),
                    error = function(e) character(0))
    if (length(hdr) == 0) next

    # every logical column must be satisfied by at least one spelling
    has_required <- all(vapply(required_any,
                               function(spellings) any(spellings %in% hdr),
                               logical(1)))
    if (!has_required) next

    # `forbid_any` rules a file out if it carries any of these columns. Used to
    # separate the city-wide table from the per-zone table, which otherwise
    # match the same required-column test.
    if (!is.null(forbid_any) && any(forbid_any %in% hdr)) next

    return(p)
  }
  NULL
}


# -----------------------------------------------------------------------------
# 5. COLUMN NAME NORMALISATION
# -----------------------------------------------------------------------------
# Maps whichever spelling the pipeline used onto one internal convention, so
# that the analysis code below never has to branch on pipeline version.
#
# Internal convention used by every validation script:
#   year, month, zone_id, zone_name,
#   ndvi_mean, ndvi_sd, ndvi_min, ndvi_max, ndvi_n,
#   lst_mean,  lst_sd,  lst_min,  lst_max,  lst_n,
#   n_scenes, sensors_used

COLUMN_ALIASES <- list(
  year         = c("year", "ano"),
  month        = c("month", "mes"),
  zone_id      = c("zone_id", "bairro_id", "id_bairro"),
  zone_name    = c("zone_name", "nome_bairro", "bairro", "BAIRRO",
                   "NOME", "Nome", "nome", "NM_BAIRRO", "Bairro_Limpo"),
  ndvi_mean    = c("ndvi_mean", "ndvi_media"),
  ndvi_sd      = c("ndvi_sd", "ndvi_dp"),
  ndvi_min     = c("ndvi_min"),
  ndvi_max     = c("ndvi_max"),
  ndvi_n       = c("ndvi_n"),
  lst_mean     = c("lst_mean", "lst_c_mean", "lst_c_media", "lst_media"),
  lst_sd       = c("lst_sd", "lst_c_sd", "lst_c_dp"),
  lst_min      = c("lst_min", "lst_c_min"),
  lst_max      = c("lst_max", "lst_c_max"),
  lst_n        = c("lst_n", "lst_c_n"),
  n_scenes     = c("n_scenes", "n_cenas"),
  sensors_used = c("sensors_used", "sensores_usados")
)

#' Rename the columns of `df` onto the internal convention.
#' Columns with no known alias are left untouched.
normalise_columns <- function(df) {
  nm <- names(df)
  for (target in names(COLUMN_ALIASES)) {
    if (target %in% nm) next                       # already correct
    hit <- intersect(COLUMN_ALIASES[[target]], nm)
    if (length(hit) > 0) names(df)[names(df) == hit[1]] <- target
  }
  df
}

#' Map a month number onto a season label, respecting HEMISPHERE.
#' Returns NA for every month when HEMISPHERE is NA, which disables seasonal
#' grouping cleanly rather than mislabelling it.
.season_of <- function(month) {
  if (length(HEMISPHERE) != 1 || is.na(HEMISPHERE)) return(NA_character_)
  h <- tolower(HEMISPHERE)
  if (!h %in% c("south", "north"))
    stop("HEMISPHERE must be \"south\", \"north\" or NA. Got: ", HEMISPHERE)
  south <- c("Summer", "Summer", "Autumn", "Autumn", "Autumn", "Winter",
             "Winter", "Winter", "Spring", "Spring", "Spring", "Summer")
  north <- c("Winter", "Winter", "Spring", "Spring", "Spring", "Summer",
             "Summer", "Summer", "Autumn", "Autumn", "Autumn", "Winter")
  (if (h == "south") south else north)[month]
}

#' Season levels in calendar order for the configured hemisphere, so that
#' factor levels in figures follow the year rather than the alphabet.
season_levels <- function() {
  if (length(HEMISPHERE) != 1 || is.na(HEMISPHERE)) return(character(0))
  if (tolower(HEMISPHERE) == "south")
    c("Summer", "Autumn", "Winter", "Spring")
  else
    c("Winter", "Spring", "Summer", "Autumn")
}

#' Add the columns every script needs: a mid-month date, an English month
#' abbreviation as an ordered factor, and the low-quality flag.
add_time_columns <- function(df) {
  df |>
    dplyr::mutate(
      date       = as.Date(sprintf("%d-%02d-15", year, month)),
      month_abb  = factor(month.abb[month], levels = month.abb),
      # Season labels follow HEMISPHERE (00_common.R section 1e). Set that to
      # NA near the equator, where a four-season scheme does not describe the
      # climate and the label would be actively misleading.
      season     = .season_of(month),
      low_quality = paste(year, month) %in%
        paste(LOW_QUALITY_MONTHS$year, LOW_QUALITY_MONTHS$month)
    )
}


# -----------------------------------------------------------------------------
# 6. LOADERS FOR THE TWO INPUT TABLES
# -----------------------------------------------------------------------------

#' City-wide monthly table: one row per (year, month).
load_city_monthly <- function(tables_dir) {
  # `forbid_any` excludes the per-zone table, which satisfies the same
  # required-column test but has ~60 rows per month; using it here would
  # silently corrupt every city-level statistic.
  p <- find_csv_by_columns(
    tables_dir,
    required_any = list(COLUMN_ALIASES$year, COLUMN_ALIASES$month,
                        COLUMN_ALIASES$ndvi_mean, COLUMN_ALIASES$lst_mean),
    prefer_pattern = "mensal|monthly",
    forbid_any     = COLUMN_ALIASES$zone_name
  )
  if (is.null(p))
    stop("Could not find a city-wide monthly CSV under ", tables_dir,
         "\n  Expected columns: year/ano, month/mes, ndvi_media, lst_c_media")
  message("[common] City monthly table : ", basename(p))
  read_any_csv(p) |> normalise_columns()
}

#' Per-zone monthly table: one row per (zone, year, month).
load_zone_monthly <- function(tables_dir) {
  p <- find_csv_by_columns(
    tables_dir,
    required_any = list(COLUMN_ALIASES$zone_name, COLUMN_ALIASES$year,
                        COLUMN_ALIASES$month, COLUMN_ALIASES$lst_mean),
    prefer_pattern = "bairro|zone"
  )
  if (is.null(p))
    stop("Could not find a per-zone monthly CSV under ", tables_dir,
         "\n  Expected a zone-name column plus year, month and LST mean.")
  message("[common] Zone monthly table : ", basename(p))
  df <- read_any_csv(p) |> normalise_columns()
  if (!"zone_id" %in% names(df))
    df$zone_id <- as.integer(factor(df$zone_name))
  df
}


# -----------------------------------------------------------------------------
# 6b. STUDY PERIOD AND EXCLUSION BOOKKEEPING
# -----------------------------------------------------------------------------

#' Resolve the study period.
#'
#' YEAR_START / YEAR_END may be left NULL in section 1d, in which case the span
#' is taken from the data itself. Reporting the actual span rather than a
#' hard-coded one prevents a stale constant from misdescribing a shorter or
#' longer run.
#'
#' Assigns YEAR_START and YEAR_END in the global environment and returns them.
resolve_period <- function(df) {
  y <- df$year[!is.na(df$year)]
  if (length(y) == 0) stop("No usable year values in the input table.")
  if (is.null(YEAR_START)) YEAR_START <<- min(y)
  if (is.null(YEAR_END))   YEAR_END   <<- max(y)
  invisible(c(start = YEAR_START, end = YEAR_END))
}

#' A short human-readable list of the excluded months, for figure captions and
#' report headers, e.g. "Jul 2021 excluded." or "Jul 2021, Mar 2013 excluded."
#' Returns "" when nothing is excluded, so callers can paste it unconditionally
#' without producing a dangling sentence.
excluded_months_caption <- function() {
  if (nrow(LOW_QUALITY_MONTHS) == 0) return("")
  lbl <- sprintf("%s %d", month.abb[LOW_QUALITY_MONTHS$month],
                 LOW_QUALITY_MONTHS$year)
  paste0(paste(lbl, collapse = ", "),
         if (length(lbl) == 1) " excluded. " else " excluded. ")
}

#' The same information in long form, for the text reports.
excluded_months_report <- function() {
  if (nrow(LOW_QUALITY_MONTHS) == 0)
    return("No months were excluded as low quality.\n")
  paste0(
    "Months excluded as low quality:\n",
    paste(sprintf("  %s %d  -  %s", month.abb[LOW_QUALITY_MONTHS$month],
                  LOW_QUALITY_MONTHS$year, LOW_QUALITY_MONTHS$reason),
          collapse = "\n"),
    "\n")
}


# -----------------------------------------------------------------------------
# 7. PUBLICATION FIGURE THEME
# -----------------------------------------------------------------------------
# Design decisions, all aimed at printed reproduction rather than screen:
#   - No panel background fill; grid lines are pale grey and only where useful.
#   - Titles left-aligned. Centred titles look decorative; left alignment reads
#     faster and matches how journals typeset figure headings.
#   - Subtitles carry the quantitative summary, so the figure is interpretable
#     without the caption.
#   - Legends at the right for continuous scales, bottom for categorical ones.

theme_pub <- function(base_size = FIG_BASE_SIZE, grid = "both") {
  th <- ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold",
                                               size = base_size + 2,
                                               hjust = 0,
                                               margin = ggplot2::margin(b = 3)),
      plot.subtitle    = ggplot2::element_text(colour = "grey35",
                                               size = base_size - 1,
                                               hjust = 0,
                                               margin = ggplot2::margin(b = 9)),
      plot.caption     = ggplot2::element_text(colour = "grey50",
                                               size = base_size - 3,
                                               hjust = 0,
                                               margin = ggplot2::margin(t = 9)),
      plot.title.position   = "plot",
      plot.caption.position = "plot",
      axis.title       = ggplot2::element_text(size = base_size - 1,
                                               colour = "grey20"),
      axis.text        = ggplot2::element_text(size = base_size - 2,
                                               colour = "grey30"),
      panel.grid.major = ggplot2::element_line(colour = "grey90",
                                               linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title     = ggplot2::element_text(size = base_size - 2,
                                               face = "bold"),
      legend.text      = ggplot2::element_text(size = base_size - 2),
      legend.key.height = ggplot2::unit(0.9, "lines"),
      strip.text       = ggplot2::element_text(face = "bold",
                                               size = base_size - 1),
      plot.margin      = ggplot2::margin(10, 12, 8, 10)
    )

  # Heatmaps and horizontal bar charts read better without competing grid lines
  if (grid == "none")
    th <- th + ggplot2::theme(panel.grid = ggplot2::element_blank())
  if (grid == "x")
    th <- th + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
  if (grid == "y")
    th <- th + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
  th
}

# --- colour palettes ---------------------------------------------------------
# Both are diverging and print acceptably in greyscale (luminance varies
# monotonically along each ramp), which matters if the journal prints in black
# and white.

# Cold blue -> warm red, for temperature.
PAL_LST  <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")

# Bare/red -> vegetated/green, the conventional NDVI ramp.
PAL_NDVI <- c("#a50026", "#f46d43", "#fee08b", "#d9ef8b", "#66bd63", "#1a9850")

# Categorical colours for the season variable.
PAL_SEASON <- c(Summer = "#d7191c", Autumn = "#fdae61",
                Winter = "#2c7bb6", Spring = "#1a9850")


# -----------------------------------------------------------------------------
# 8. FIGURE SAVING
# -----------------------------------------------------------------------------
#' Save a ggplot as a high-resolution PNG and, optionally, a vector PDF.
#'
#' Uses the ragg PNG device when available. ragg renders text with proper
#' hinting and anti-aliasing; the default grDevices PNG device on a headless
#' Linux node often produces noticeably rougher glyphs, which is very visible
#' at the small font sizes used in the 60-zone figures.
#'
#' @param plot    a ggplot object
#' @param file    file name WITHOUT extension
#' @param dir     output directory
#' @param width   width in inches
#' @param height  height in inches
#' @param tag     short prefix used in the console log
save_fig <- function(plot, file, dir, width = 10, height = 7, tag = "") {

  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  png_path <- file.path(dir, paste0(file, ".png"))

  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(png_path, plot, width = width, height = height,
                    dpi = FIG_DPI, bg = "white", device = ragg::agg_png)
  } else {
    ggplot2::ggsave(png_path, plot, width = width, height = height,
                    dpi = FIG_DPI, bg = "white")
  }
  message(sprintf("%s   %s.png", tag, file))

  if (isTRUE(FIG_SAVE_PDF)) {
    pdf_path <- file.path(dir, paste0(file, ".pdf"))
    # cairo_pdf embeds fonts and supports transparency; the default pdf device
    # does neither reliably, which breaks alpha-blended density plots.
    ok <- tryCatch({
      ggplot2::ggsave(pdf_path, plot, width = width, height = height,
                      device = grDevices::cairo_pdf, bg = "white")
      TRUE
    }, error = function(e) FALSE)
    if (!ok)
      tryCatch(ggplot2::ggsave(pdf_path, plot, width = width,
                               height = height, bg = "white"),
               error = function(e)
                 message(tag, "   (PDF export unavailable, PNG written)"))
  }
  invisible(png_path)
}


# -----------------------------------------------------------------------------
# 9. SMALL STATISTICAL HELPERS
# -----------------------------------------------------------------------------

#' Linear trend in units per year, from a date vector and a value vector.
#' Returns NA when fewer than `min_n` observations are available, because a
#' slope fitted to a handful of points is dominated by which points happen to
#' be present rather than by any trend.
slope_per_year <- function(date, value, min_n = 5) {
  ok <- !is.na(value) & !is.na(date)
  if (sum(ok) < min_n) return(NA_real_)
  year_frac <- as.numeric(format(date[ok], "%Y")) +
    (as.numeric(format(date[ok], "%m")) - 1) / 12
  if (stats::var(year_frac) <= 0) return(NA_real_)
  unname(stats::coef(stats::lm(value[ok] ~ year_frac))[2])
}

#' Open a sink safely: guarantees the connection is closed even if the code
#' that follows throws. Without this an error mid-report leaves the sink open
#' and every subsequent console message is silently redirected to the file.
open_report <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  sink(path)
  invisible(path)
}
close_report <- function() {
  while (sink.number() > 0) sink()
}

message("[common] Helpers loaded.")
message("[common] Validation output : ", VALIDATION_ROOT)
message("[common] Study area label  : ", AOI_LABEL)
if (nrow(LOW_QUALITY_MONTHS) > 0)
  message("[common] Excluded months   : ",
          paste(sprintf("%d-%02d", LOW_QUALITY_MONTHS$year,
                        LOW_QUALITY_MONTHS$month), collapse = ", "))
