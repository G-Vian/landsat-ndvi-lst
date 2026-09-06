# =============================================================================
# MODULE 04 - SPATIAL CLIPPING TO THE AREA OF INTEREST (AOI)
# =============================================================================
# (This module was formerly 04_clip_santos.R. It is now area-agnostic: the AOI
#  comes entirely from SHP_FILE in 00_config.R.)
#
# Responsibility:
#   - Load the AOI polygon file and derive two objects from it:
#       * the OUTER BOUNDARY (all polygons dissolved) -> used to clip rasters
#       * the ZONES (each polygon kept separate)      -> used for zonal stats
#   - Reproject both AOI and rasters to CRS_TARGET
#   - Clip (crop + mask) NDVI/LST rasters to the AOI
#   - Save the clipped rasters as GeoTIFF
#
# YOUR SHAPEFILE IS THE ONLY THING THAT DEFINES THE STUDY AREA.
# Requirements are documented in 00_config.R Section 4. In short: valid CRS,
# non-overlapping polygons, and a name column for the zones.
# =============================================================================

source("00_config.R")

# Session cache so the shapefile is read from disk only once
.aoi_cache   <- NULL
.zones_cache <- NULL

# -----------------------------------------------------------------------------
# 1. LOADING THE AOI
# -----------------------------------------------------------------------------

#' Load the AOI outer boundary: all polygons dissolved into a single geometry.
#'
#' Used for clipping, not for statistics. Dissolving matters when the input has
#' many polygons (neighbourhoods, municipalities): we want one outline, and we
#' want it valid, hence st_make_valid() after the union.
#'
#' @return sf object with a single (multi)polygon in CRS_TARGET
load_aoi_boundary <- function() {

  if (!is.null(.aoi_cache)) return(.aoi_cache)

  if (!file.exists(SHP_FILE))
    stop("[04_clip] AOI file not found: ", SHP_FILE,
         "\n  Set SHP_FILE in 00_config.R Section 4.")

  shp <- sf::st_read(SHP_FILE, quiet = TRUE)

  if (is.na(sf::st_crs(shp)))
    stop("[04_clip] The AOI file has NO CRS defined. ",
         "Assign one in GIS software before using it — the pipeline can ",
         "reproject, but it cannot guess the original projection.")

  log_msg(sprintf("[04_clip] AOI loaded: %d polygon(s) | CRS: %s",
                  nrow(shp), sf::st_crs(shp)$input))

  boundary <- shp |>
    sf::st_union() |>
    sf::st_as_sf() |>
    sf::st_make_valid() |>
    sf::st_transform(CRS_TARGET)

  log_msg(sprintf("[04_clip] Boundary dissolved and reprojected to %s.", CRS_TARGET))

  .aoi_cache <<- boundary
  boundary
}


#' Load the individual AOI zones (one row per polygon) for zonal statistics.
#'
#' Adds a `zone_name` column resolved from ZONE_NAME_COLUMN (00_config.R), or
#' auto-detected when that is NULL.
#'
#' @return sf object in CRS_TARGET, with a guaranteed `zone_name` column
load_aoi_zones <- function() {

  if (!is.null(.zones_cache)) return(.zones_cache)

  if (!file.exists(SHP_FILE))
    stop("[04_clip] AOI file not found: ", SHP_FILE)

  shp <- sf::st_read(SHP_FILE, quiet = TRUE) |>
    sf::st_make_valid() |>
    sf::st_transform(CRS_TARGET)

  name_col <- resolve_zone_name_column(shp)
  shp$zone_name <- as.character(sf::st_drop_geometry(shp)[[name_col]])

  # Duplicate names would silently merge distinct zones in later group_by()
  dups <- shp$zone_name[duplicated(shp$zone_name)]
  if (length(dups) > 0)
    log_msg(sprintf(
      "[04_clip] WARNING: duplicated zone names in column '%s': %s. ",
      name_col, paste(unique(dups), collapse = ", ")), "WARN")

  log_msg(sprintf("[04_clip] %d zones loaded (name column: '%s').",
                  nrow(shp), name_col))

  .zones_cache <<- shp
  shp
}


#' Decide which attribute column holds the zone name.
#'
#' Explicit configuration wins. Auto-detection is a convenience for quick
#' tests, but it can pick the wrong column on a rich shapefile — set
#' ZONE_NAME_COLUMN in 00_config.R for anything you intend to publish.
resolve_zone_name_column <- function(shp) {

  cols <- colnames(sf::st_drop_geometry(shp))

  # (1) Explicit configuration
  if (!is.null(ZONE_NAME_COLUMN)) {
    if (!ZONE_NAME_COLUMN %in% cols)
      stop(sprintf(
        "[04_clip] ZONE_NAME_COLUMN = '%s' is not a column of the AOI file.\n  Available columns: %s",
        ZONE_NAME_COLUMN, paste(cols, collapse = ", ")))
    return(ZONE_NAME_COLUMN)
  }

  # (2) Common conventions: generic English, Brazilian IBGE, GADM
  candidates <- c("NOME", "Nome", "nome",
                  "NAME", "Name", "name",
                  "NM_BAIRRO", "nm_bairro", "BAIRRO", "bairro",
                  "NM_MUN", "NM_MUNICIP", "MUNICIPIO", "municipio",
                  "NAME_1", "NAME_2", "NAME_3",          # GADM levels
                  "district", "DISTRICT", "zone", "ZONE")
  found <- candidates[candidates %in% cols]
  if (length(found) > 0) {
    log_msg(sprintf("[04_clip] Zone name column auto-detected: '%s'. ",
                    found[1]), "WARN")
    return(found[1])
  }

  # (3) Fallback: first character column
  chr_cols <- cols[sapply(sf::st_drop_geometry(shp), is.character)]
  if (length(chr_cols) > 0) {
    log_msg(sprintf(
      "[04_clip] No standard name column found; falling back to '%s'. ",
      chr_cols[1]), "WARN")
    return(chr_cols[1])
  }

  stop("[04_clip] Could not determine a zone name column. ",
       "Set ZONE_NAME_COLUMN explicitly in 00_config.R.")
}

# -----------------------------------------------------------------------------
# 2. CLIPPING
# -----------------------------------------------------------------------------

#' Reproject and clip a raster to the AOI.
#'
#' Order of operations matters: NDVI and LST are computed FIRST, in the native
#' projection, and only then reprojected. Reprojecting the raw bands before the
#' index would interpolate DNs across cloud edges and contaminate the result.
#'
#' Bilinear interpolation is used because NDVI and LST are continuous. Use
#' method = "near" only if you ever adapt this to a categorical product.
#'
#' @param r        SpatRaster to clip
#' @param boundary sf AOI outline (defaults to load_aoi_boundary())
#' @return clipped and masked SpatRaster
clip_to_aoi <- function(r, boundary = load_aoi_boundary()) {

  crs_r      <- terra::crs(r, describe = TRUE)$code
  crs_target <- gsub("EPSG:", "", CRS_TARGET)

  if (!isTRUE(crs_r == crs_target))
    r <- terra::project(r, CRS_TARGET, method = "bilinear")

  boundary_v <- terra::vect(boundary)

  # crop() reduces to the bounding box (fast); mask() applies the actual
  # polygon shape (sets pixels outside the geometry to NA).
  r_crop <- terra::crop(r, boundary_v)
  terra::mask(r_crop, boundary_v)
}

# -----------------------------------------------------------------------------
# 3. SAVING CLIPPED RASTERS
# -----------------------------------------------------------------------------

#' Write clipped NDVI and LST rasters as compressed GeoTIFFs.
#'
#' FLT4S (32-bit float) preserves the decimal precision of both products.
#' LZW compression is lossless. Expect roughly 1-5 MB per scene for a
#' city-sized AOI; scale that estimate up linearly with AOI area and set
#' SAVE_CLIPPED_RASTERS = FALSE in 00_config.R if disk becomes the constraint.
#'
#' @return list(ndvi = <path>, lst = <path>)
save_rasters <- function(ndvi_r, lst_r, scene_id) {

  dir_ndvi <- file.path(OUTPUT_DIR, "rasters", "ndvi")
  dir_lst  <- file.path(OUTPUT_DIR, "rasters", "lst")
  dir.create(dir_ndvi, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_lst,  recursive = TRUE, showWarnings = FALSE)

  out_ndvi <- file.path(dir_ndvi, paste0(scene_id, "_NDVI.tif"))
  out_lst  <- file.path(dir_lst,  paste0(scene_id, "_LST_Celsius.tif"))

  terra::writeRaster(ndvi_r, out_ndvi, overwrite = TRUE,
                     datatype = "FLT4S",
                     gdal = c("COMPRESS=LZW", "TILED=YES"))
  terra::writeRaster(lst_r, out_lst, overwrite = TRUE,
                     datatype = "FLT4S",
                     gdal = c("COMPRESS=LZW", "TILED=YES"))

  log_msg(sprintf("[04_clip] Rasters saved: %s", scene_id))
  list(ndvi = out_ndvi, lst = out_lst)
}

# -----------------------------------------------------------------------------
# 4. FOOTPRINT OVERLAP CHECK
# -----------------------------------------------------------------------------

#' Does this scene cover enough of the AOI to be worth processing?
#'
#' ############################## IMPORTANT ##############################
#' The threshold MIN_AOI_OVERLAP_FRAC (00_config.R Section 8c) is the single
#' setting most likely to break when you move to a LARGER study area.
#'
#' One Landsat scene is ~185 x 180 km. For a city, one scene covers 100% of
#' the AOI and the default 0.10 works. For a state or country, each scene
#' covers only a small slice, so a 0.10 threshold rejects EVERYTHING and the
#' run finishes with zero scenes.
#'
#' If your run reports "no overlap" for every scene, this is why.
#' #######################################################################
#'
#' @return TRUE if the intersection exceeds MIN_AOI_OVERLAP_FRAC of the AOI
has_overlap <- function(r, boundary = load_aoi_boundary()) {

  # Threshold of 0 (or less) disables the check entirely
  if (MIN_AOI_OVERLAP_FRAC <= 0) return(TRUE)

  # Guard against rasters with an empty CRS, which would crash st_set_crs
  crs_r_str <- terra::crs(r)
  if (is.na(crs_r_str) || nchar(trimws(crs_r_str)) == 0)
    crs_r_str <- CRS_TARGET

  ext_r <- sf::st_as_sf(terra::as.polygons(terra::ext(r))) |>
    sf::st_set_crs(crs_r_str) |>
    sf::st_transform(sf::st_crs(boundary))

  inter <- suppressWarnings(sf::st_intersection(boundary, ext_r))

  if (nrow(inter) == 0 || all(sf::st_is_empty(inter))) return(FALSE)

  coverage <- as.numeric(sum(sf::st_area(inter)) / sum(sf::st_area(boundary)))
  coverage > MIN_AOI_OVERLAP_FRAC
}

# -----------------------------------------------------------------------------
# BACKWARD-COMPATIBILITY ALIASES
# -----------------------------------------------------------------------------
# Kept so that scripts written against the Santos-specific version still run.
load_santos_boundary <- load_aoi_boundary
load_santos_bairros  <- load_aoi_zones
clip_to_santos       <- clip_to_aoi
