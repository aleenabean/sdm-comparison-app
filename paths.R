# centralized path definitions for project

PROJECT_ROOT <- normalizePath(".", winslash = "/", mustWork = FALSE)

# top-level directories
DIR_RAW          <- file.path(PROJECT_ROOT, "rdata_raw")
DIR_INTERMEDIATE <- file.path(PROJECT_ROOT, "data_intermediate")
DIR_PROCESSED    <- file.path(PROJECT_ROOT, "data_processed")
DIR_MODELS       <- file.path(PROJECT_ROOT, "models")
DIR_EVAL         <- file.path(PROJECT_ROOT, "eval")

# subdirectories in rdata_raw
DIR_RAW_INAT     <- file.path(DIR_RAW, "inat")
DIR_RAW_IUCN     <- file.path(DIR_RAW, "iucn")
DIR_RAW_WORLDCLIM<- file.path(DIR_RAW, "worldclim")


# function to check that required directories exist
check_project_dirs <- function() {
  dirs <- c(
    DIR_RAW, DIR_INTERMEDIATE, DIR_PROCESSED,
    DIR_MODELS, DIR_EVAL,
    DIR_RAW_INAT, DIR_RAW_IUCN, DIR_RAW_WORLDCLIM
  )
  missing <- dirs[!dir.exists(dirs)]
  if (length(missing) > 0) {
    warning("These directories are missing:\n",
            paste(missing, collapse = "\n"),
            "\nRun R/00_setup_project.R to create them.")
  } else {
    message("All expected directories exist.")
  }
}


