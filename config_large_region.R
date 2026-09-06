# =============================================================================
# EXAMPLE CONFIG - LARGE REGION (STATE / PROVINCE / COUNTRY)
# =============================================================================
# Profile: an AOI spanning MANY Landsat scenes, e.g. a whole state.
# This differs most from a city-scale run, and is where the defaults will
# actively break if left unchanged.
#
# THE FOUR THINGS THAT MUST CHANGE:
#   1. MIN_AOI_OVERLAP_FRAC  -> near zero (each scene covers a small slice)
#   2. COMPUTE_PIXEL_TABLES  -> FALSE (hundreds of millions of rows otherwise)
#   3. CRS_TARGET            -> equal-area, not a single UTM zone
#   4. SAVE_CLIPPED_RASTERS  -> consider FALSE if disk is limited
# =============================================================================

AOI_NAME   <- "MyState"
BASE_DIR   <- "/scratch/user/state_project"
ZONE_LABEL <- "municipality"
ZONE_NAME_COLUMN <- "NM_MUN"

DIRS_LANDSAT <- list(all = file.path(BASE_DIR, "landsat"))
SHP_FILE     <- file.path(BASE_DIR, "shp/municipalities.shp")
OUTPUT_DIR   <- file.path(BASE_DIR, "results")

# A large AOI multiplies the scene count enormously (many path/rows x dates).
# Start with a SHORT period to size the run before committing to the full one.
YEAR_START <- 2020
YEAR_END   <- 2024
TIER_KEEP  <- "T1"

# (1) EQUAL-AREA CRS: a region spanning several UTM zones is distorted by any
# single UTM projection. Since zonal statistics weight by pixel count, and
# pixel ground area varies with that distortion, this would bias zone
# comparisons.
#   South America Albers .. "ESRI:102033"
#   CONUS Albers .......... "EPSG:5070"
#   Europe LAEA ........... "EPSG:3035"
CRS_TARGET <- "ESRI:102033"

# A large region spans several climates: widen the bounds so no sub-region is
# systematically truncated, then narrow them after inspecting
# tables/anomalies/ from the first run.
LST_MIN <- -10.0
LST_MAX <-  75.0

QA_BITS_MASK <- c(0, 1, 3, 4, 5, 7)

# (2) THE CRITICAL SETTING. One scene covers a tiny fraction of a state, so
# the default 0.10 would reject EVERY scene and the run would end empty.
MIN_AOI_OVERLAP_FRAC <- 0.0001      # or 0 to disable the check entirely

# The AOI-wide validity check is also near-meaningless at this scale: a single
# scene never covers most of the region. Keep it very low.
MIN_AOI_VALID_FRAC   <- 0.001
MIN_SCENE_VALID_FRAC <- 0.10        # scene-level cloud filter still useful

# (3) SCALING: per-pixel tables are NOT tractable at this size.
SAVE_CLIPPED_RASTERS <- TRUE        # set FALSE if disk-constrained
COMPUTE_ZONE_STATS   <- TRUE        # scales with polygon count, safe
COMPUTE_PIXEL_TABLES <- FALSE       # MUST stay FALSE
GENERATE_PLOTS       <- TRUE        # consider FALSE: one PNG per period

# (4) Lower terra's memory ceiling on a shared cluster node
TERRA_MEMFRAC <- 0.4

# NOTE ON MOSAICKING: this pipeline processes each scene independently and
# aggregates statistically (weighted by valid pixel count). It does NOT build
# a seamless mosaic. For a large AOI a given zone may be covered by different
# scenes on different dates - handled correctly for zonal statistics, but the
# per-pixel maps of module 08 can show scene-boundary seams.
