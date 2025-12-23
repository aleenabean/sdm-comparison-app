# helper functions to prepare presence + pseudo-absence data
# for a single-species GLM SDM using iNat + WorldClim

library(dplyr)
library(readr)
library(terra)
library(googledrive)

source("paths.R")
check_project_dirs()

FILE_ENV_RDS  <- file.path(DIR_INTERMEDIATE, "worldclim21_bio_2-5m_stack.rds")
FILE_INAT_CSV <- file.path(DIR_INTERMEDIATE, "inat_occurrences_1000_species.csv")

# Google Drive file ID for the zipped WorldClim folder
# This zip contains the WorldClim .tif files (ex wc2.1_2.5m_bio_1.tif, ...)
WORLDCLIM_DRIVE_ID <- "1iYwuMdbGTSMgPuY0VYMPBZq6kP2TPSLk"


# WorldClim download + env stack loader
ensure_worldclim_local <- function() {
  existing_tifs <- list.files(
    DIR_RAW_WORLDCLIM,
    pattern   = "\\.tif$",
    full.names = TRUE,
    recursive  = TRUE
  )
  
  if (length(existing_tifs) > 0) {
    return(invisible(NULL))
  }
  
  # otherwise get from Google Drive
  dir.create(DIR_RAW_WORLDCLIM, recursive = TRUE, showWarnings = FALSE)
  
  #download zip to a temp file
  tmp_zip <- file.path(tempdir(), "worldclim_wc2_1_2_5m_bio.zip")
  
  cat("WorldClim .tif files not found; downloading zip from Google Drive...\n")
  googledrive::drive_deauth()
  googledrive::drive_download(
    googledrive::as_id(WORLDCLIM_DRIVE_ID),
    path      = tmp_zip,
    overwrite = TRUE
  )
  
  cat("Unzipping WorldClim data...\n")
  tmp_unzip_dir <- file.path(tempdir(), "worldclim_unzip")
  dir.create(tmp_unzip_dir, recursive = TRUE, showWarnings = FALSE)
  utils::unzip(tmp_zip, exdir = tmp_unzip_dir)
  
  # find all .tif files in the unzipped content and copy them into DIR_RAW_WORLDCLIM
  unzipped_tifs <- list.files(
    tmp_unzip_dir,
    pattern   = "\\.tif$",
    full.names = TRUE,
    recursive  = TRUE
  )
  
  if (length(unzipped_tifs) == 0) {
    stop(
      "No .tif files found inside the downloaded WorldClim zip.\n",
      "Check that the Google Drive file (ID = ", WORLDCLIM_DRIVE_ID, ") ",
      "contains the WorldClim .tif rasters."
    )
  }
  
  file.copy(unzipped_tifs, DIR_RAW_WORLDCLIM, overwrite = TRUE)
  
  # clean up temp files
  unlink(tmp_zip)
  unlink(tmp_unzip_dir, recursive = TRUE)
  
  cat("WorldClim download + unzip complete. Stored .tif files in:\n  ",
      DIR_RAW_WORLDCLIM, "\n")
}

# load the WorldClim environmental stack from local .tif files
load_env_stack_default <- function() {
  ensure_worldclim_local()
  
  tif_files <- list.files(
    DIR_RAW_WORLDCLIM,
    pattern   = "\\.tif$",
    full.names = TRUE,
    recursive  = TRUE
  )
  
  if (length(tif_files) == 0) {
    stop(
      "No WorldClim .tif files found in DIR_RAW_WORLDCLIM: ", DIR_RAW_WORLDCLIM,
      "\nDownload failed or data is missing."
    )
  }
  
  terra::rast(tif_files)
}

load_inat_occurrences_default <- function(path = FILE_INAT_CSV) {
  if (!file.exists(path)) {
    stop("inat_occurrences_1000_species.csv not found at: ", path,
         "\nRun the iNat download/prep scripts or adjust FILE_INAT_CSV.")
  }
  readr::read_csv(path, show_col_types = FALSE)
}

# helper to guess column names

guess_species_column <- function(df) {
  candidates <- c(
    "iucn_scientific_name",  
    "scientific_name",
    "sci_name",
    "species",
    "species_name",
    "scientificName"         # GBIF column if nothing else exists
  )
  
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) {
    stop("Could not find a species-name column in inat_occurrences.\n",
         "Looked for: ", paste(candidates, collapse = ", "), "\n",
         "Columns are: ", paste(names(df), collapse = ", "))
  }
  hit[1]
}

guess_lon_lat_columns <- function(df) {
  lon_candidates <- c("lon", "longitude", "decimalLongitude")
  lat_candidates <- c("lat", "latitude", "decimalLatitude")
  
  lon_hit <- intersect(lon_candidates, names(df))
  lat_hit <- intersect(lat_candidates, names(df))
  
  if (length(lon_hit) == 0 || length(lat_hit) == 0) {
    stop("Could not find lon/lat columns in inat_occurrences.\n",
         "Looked for lon: ", paste(lon_candidates, collapse = ", "),
         " and lat: ", paste(lat_candidates, collapse = ", "), "\n",
         "Columns are: ", paste(names(df), collapse = ", "))
  }
  
  list(lon = lon_hit[1], lat = lat_hit[1])
}

deg2rad <- function(x) x * pi / 180


# main: prepare presence + pseudo-absence data for one species
# returns a data.frame with columns: 
#presence (0/1), lon, lat, x1..x4, and all WorldClim variables

# arguments:
# scientific_name: character, species name to filter
# n_pseudo_abs: max number of pseudo-absence points to sample
# env_stack: optional SpatRaster; if NULL, loads from Google Drive/WorldClim
# inat_occ: optional data.frame; if NULL, loads FILE_INAT_CSV
# bbox_margin: degrees to expand presence bounding box

prepare_glm_data_for_species <- function(
    scientific_name,
    n_pseudo_abs = 1000L,
    env_stack    = NULL,
    inat_occ     = NULL,
    bbox_margin  = 1
) {
  
  if (is.null(env_stack)) {
    env_stack <- load_env_stack_default()
  }
  if (is.null(inat_occ)) {
    inat_occ <- load_inat_occurrences_default()
  }
  
  sp_col  <- guess_species_column(inat_occ)
  ll_cols <- guess_lon_lat_columns(inat_occ)
  lon_col <- ll_cols$lon
  lat_col <- ll_cols$lat
  
  pres_raw <- inat_occ %>%
    filter(.data[[sp_col]] == scientific_name)
  
  if (nrow(pres_raw) == 0) {
    stop("No iNat/GBIF occurrences found for species: ", scientific_name)
  }
  
  pres_raw <- pres_raw %>%
    mutate(
      lon = as.numeric(.data[[lon_col]]),
      lat = as.numeric(.data[[lat_col]])
    ) %>%
    filter(!is.na(lon), !is.na(lat))
  
  if (nrow(pres_raw) == 0) {
    stop("All lon/lat NA after coercion for species: ", scientific_name)
  }
  
  lon_min <- min(pres_raw$lon) - bbox_margin
  lon_max <- max(pres_raw$lon) + bbox_margin
  lat_min <- max(min(pres_raw$lat) - bbox_margin, -90)
  lat_max <- min(max(pres_raw$lat) + bbox_margin,  90)
  
  lon_min <- max(lon_min, -180)
  lon_max <- min(lon_max,  180)
  
  pts_pres <- terra::vect(
    pres_raw %>% select(lon, lat),
    geom = c("lon", "lat"),
    crs  = "EPSG:4326"
  )
  
  env_pres <- terra::extract(env_stack, pts_pres, df = TRUE)
  env_pres <- env_pres[, -1, drop = FALSE]  # drop ID column
  env_pres[is.na(env_pres)] <- 0
  
  lon_rad_pres <- deg2rad(pres_raw$lon)
  lat_rad_pres <- deg2rad(pres_raw$lat)
  
  pres_df <- data.frame(
    presence = 1L,
    lon      = pres_raw$lon,
    lat      = pres_raw$lat,
    x1       = sin(lon_rad_pres),
    x2       = cos(lon_rad_pres),
    x3       = sin(lat_rad_pres),
    x4       = cos(lat_rad_pres),
    env_pres
  )
  
  # pseudo-absences (background)
  n_pres <- nrow(pres_df)
  
  #at most 10x presences
  n_bg_target <- min(n_pseudo_abs, n_pres * 10L)
  
  if (n_bg_target < 1L) {
    warning("n_bg_target < 1; returning presences only.")
    return(pres_df)
  }
  
  set.seed(123)  
  lon_bg <- runif(n_bg_target, min = lon_min, max = lon_max)
  lat_bg <- runif(n_bg_target, min = lat_min, max = lat_max)
  
  bg_raw <- data.frame(lon = lon_bg, lat = lat_bg)
  
  pts_bg <- terra::vect(
    bg_raw,
    geom = c("lon", "lat"),
    crs  = "EPSG:4326"
  )
  
  env_bg <- terra::extract(env_stack, pts_bg, df = TRUE)
  env_bg <- env_bg[, -1, drop = FALSE]
  env_bg[is.na(env_bg)] <- 0
  
  lon_rad_bg <- deg2rad(bg_raw$lon)
  lat_rad_bg <- deg2rad(bg_raw$lat)
  
  bg_df <- data.frame(
    presence = 0L,
    lon      = bg_raw$lon,
    lat      = bg_raw$lat,
    x1       = sin(lon_rad_bg),
    x2       = cos(lon_rad_bg),
    x3       = sin(lat_rad_bg),
    x4       = cos(lat_rad_bg),
    env_bg
  )
  
  # combine & return
  glm_df <- dplyr::bind_rows(pres_df, bg_df)
  
  glm_df <- glm_df %>%
    mutate(
      presence = as.integer(presence),
      lon      = as.numeric(lon),
      lat      = as.numeric(lat)
    )
  
  message("Prepared GLM data for ", scientific_name, ": ",
          n_pres, " presences + ", nrow(bg_df), " pseudo-absences (",
          nrow(glm_df), " rows total).")
  
  glm_df
}
