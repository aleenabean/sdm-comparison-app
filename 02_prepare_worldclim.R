#read WorldClim 2.1 bioclim (2.5 arc-min)

library(terra)

source("paths.R")
check_project_dirs()


bio_dirs <- list.dirs(DIR_RAW_WORLDCLIM, full.names = TRUE, recursive = FALSE)
bio_dirs <- bio_dirs[grepl("bio", basename(bio_dirs), ignore.case = TRUE)]

bio_zips <- list.files(
  DIR_RAW_WORLDCLIM,
  pattern   = "\\.zip$",
  full.names = TRUE
)

if (length(bio_dirs) == 0 && length(bio_zips) == 0) {
  stop(
    "Could not find WorldClim bioclim data in ",
    DIR_RAW_WORLDCLIM,
    "\nExpected either a folder like 'wc2.1_2.5m_bio' or a zip file."
  )
}

#since we're switching around with files and zipes
if (length(bio_dirs) > 0) {
  bio_source_type <- "dir"
  bio_source_path <- bio_dirs[1]
  message("Using WorldClim bioclim directory: ", bio_source_path)
  
  tif_dir <- bio_source_path
} else {
  bio_source_type <- "zip"
  bio_source_path <- bio_zips[1]
  message("Using WorldClim bioclim zip: ", bio_source_path)
  
  tmp_dir <- tempfile(pattern = "worldclim_bio_")
  dir.create(tmp_dir)
  unzip(bio_source_path, exdir = tmp_dir)
  tif_dir <- tmp_dir
}

#loading tif files

tif_files <- list.files(tif_dir, pattern = "\\.tif$", full.names = TRUE)
if (length(tif_files) == 0) {
  stop("No .tif files found in: ", tif_dir)
}

message("Found ", length(tif_files), " WorldClim bioclim tif(s).")

bio_stack <- rast(tif_files)
message("WorldClim bioclim stack has ", nlyr(bio_stack), " layers.")
message("Layer names:")
print(names(bio_stack))


env_stack <- bio_stack

#save rds for later

file_env_rds <- file.path(DIR_INTERMEDIATE, "worldclim21_bio_2-5m_stack.rds")
saveRDS(env_stack, file_env_rds)

message("Saved WorldClim env stack to: ", file_env_rds)
