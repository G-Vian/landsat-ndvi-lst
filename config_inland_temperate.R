# =============================================================================
# EXAMPLE CONFIG - INLAND TEMPERATE CITY WITH WINTER FROST
# =============================================================================
# Profile: mid-latitude continental city with sub-zero winters and no coast.
# Illustrates the settings that MUST change for a colder, drier climate.
# =============================================================================

AOI_NAME   <- "InlandCity"
BASE_DIR   <- "/home/user/research/inland_city"
ZONE_LABEL <- "district"
ZONE_NAME_COLUMN <- "district_name"

DIRS_LANDSAT <- list(all = file.path(BASE_DIR, "landsat"))
SHP_FILE     <- file.path(BASE_DIR, "shp/districts.gpkg")
OUTPUT_DIR   <- file.path(BASE_DIR, "results")

YEAR_START <- 2013     # Landsat 8 onwards: one consistent sensor generation
YEAR_END   <- 2025
TIER_KEEP  <- "T1"

# WGS 84 / UTM zone 33N - replace with the zone covering your city
CRS_TARGET <- "EPSG:32633"

# CRITICAL CHANGE: frost occurs here. The coastal default of 5 C would
# discard every genuine winter observation.
LST_MIN <- -25.0
LST_MAX <-  60.0

# Snow bit matters far more here than at tropical latitudes.
# Water bit kept: rivers and lakes would still distort land statistics.
QA_BITS_MASK <- c(0, 1, 3, 4, 5, 7)

# Winter scenes are frequently snow- or cloud-covered; relaxing the
# scene-level threshold keeps partially usable acquisitions.
MIN_SCENE_VALID_FRAC <- 0.05
MIN_AOI_VALID_FRAC   <- 0.15
MIN_AOI_OVERLAP_FRAC <- 0.10

SAVE_CLIPPED_RASTERS <- TRUE
COMPUTE_ZONE_STATS   <- TRUE
COMPUTE_PIXEL_TABLES <- TRUE
GENERATE_PLOTS       <- TRUE

# Also edit in R/08_plots.R so the colour scale spans the winter range
# (or leave both NULL to let the pipeline choose them from your data):
#   PLOT_LIMITS_LST <- c(-10, 45)
