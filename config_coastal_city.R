# =============================================================================
# EXAMPLE CONFIG - SMALL COASTAL CITY  (the original Santos-SP setup)
# =============================================================================
# Profile: humid subtropical coastal municipality, ~40 km2, tens of
# neighbourhoods. One Landsat scene fully covers the AOI.
# Copy the blocks below into R/00_config.R.
# =============================================================================

AOI_NAME   <- "Santos-SP"
BASE_DIR   <- "/home/user/research/santos"
ZONE_LABEL <- "neighbourhood"
ZONE_NAME_COLUMN <- "BAIRRO"

DIRS_LANDSAT <- list(
  l04_l05 = file.path(BASE_DIR, "landsat_data/l04_l05"),
  l07     = file.path(BASE_DIR, "landsat_data/l07"),
  l08_l09 = file.path(BASE_DIR, "landsat_data/l08_l09")
)
SHP_FILE   <- file.path(BASE_DIR, "shp/neighbourhoods.shp")
OUTPUT_DIR <- file.path(BASE_DIR, "results_landsat")

YEAR_START  <- 2010
YEAR_END    <- 2025
MONTHS_KEEP <- NULL
TIER_KEEP   <- "T1"

# SIRGAS 2000 / UTM 23S - correct for the Sao Paulo coast
CRS_TARGET <- "EPSG:31983"

# Coastal subtropical: never freezes; dark roofs reach ~70 C in summer.
# These are the values used in the original study, with the reasoning
# recorded there: 5 C = coldest plausible coastal winter surface, below
# which pixels are shadow or scene border; 70 C = metal roofing under
# extreme summer insolation, above which values are artefactual.
LST_MIN <- 5.0
LST_MAX <- 70.0

# Water masked: the AOI is an island, and sea-surface LST would swamp the
# urban-heat signal in coastal neighbourhoods.
QA_BITS_MASK <- c(0, 1, 3, 4, 5, 7)

# AOI smaller than a single scene -> the default overlap rule works
MIN_AOI_OVERLAP_FRAC <- 0.10
MIN_SCENE_VALID_FRAC <- 0.10
MIN_AOI_VALID_FRAC   <- 0.10

# Small AOI -> per-pixel tables are tractable (~44k rows/month)
SAVE_CLIPPED_RASTERS <- TRUE
COMPUTE_ZONE_STATS   <- TRUE
COMPUTE_PIXEL_TABLES <- TRUE
GENERATE_PLOTS       <- TRUE
