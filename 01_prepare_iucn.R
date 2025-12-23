# Read IUCN mammals + amphibians shapefiles, combine them, 
#and build a species table

library(sf)
library(dplyr)
library(readr)

source("paths.R")
check_project_dirs()

read_iucn_zip <- function(zip_path) {
  if (!file.exists(zip_path)) {
    stop("Zip file not found: ", zip_path)
  }
  
# unzip
  tmp_dir <- tempfile(pattern = "iucn_unzip_")
  dir.create(tmp_dir)
  unzip(zip_path, exdir = tmp_dir)
  
# list all shapefiles
  shp_files <- list.files(tmp_dir, pattern = "\\.shp$", full.names = TRUE)
  if (length(shp_files) == 0) {
    stop("No .shp files found inside: ", zip_path)
  }
  
  message("Found ", length(shp_files), " shapefile(s) in ", zip_path)
  
#read and row-bind them
  sf_list <- lapply(shp_files, function(f) {
    message("Reading shapefile: ", f)
    st_read(f, quiet = TRUE)
  })
  
  do.call(dplyr::bind_rows, sf_list)
}


zip_mammals    <- file.path(DIR_RAW_IUCN, "mammals_ranges.zip")
zip_amphibians <- file.path(DIR_RAW_IUCN, "amphibians_ranges.zip")

iucn_mammals    <- read_iucn_zip(zip_mammals)
iucn_amphibians <- read_iucn_zip(zip_amphibians)

#combine mammals + amphibians into one sf object

iucn_all <- bind_rows(
  iucn_mammals %>% mutate(taxon_group = "mammals"),
  iucn_amphibians %>% mutate(taxon_group = "amphibians")
)

message("Columns in IUCN data:")
print(names(iucn_all))

name_col <- "sci_name"
message("Using species name column: ", name_col)

iucn_species_table <- iucn_all %>%
  st_drop_geometry() %>%
  select(taxon_group, scientific_name = !!sym(name_col)) %>%
  distinct() %>%
  arrange(taxon_group, scientific_name) %>%
  mutate(
    species_index = row_number() - 1L 
  )


#save data
file_ranges_rds <- file.path(DIR_INTERMEDIATE, "iucn_ranges_mammals_amphibians.rds")
file_species_csv <- file.path(DIR_INTERMEDIATE, "iucn_species_table_mammals_amphibians.csv")

saveRDS(iucn_all, file_ranges_rds)
write_csv(iucn_species_table, file_species_csv)

message("Saved combined IUCN ranges to: ", file_ranges_rds)
message("Saved IUCN species table to: ", file_species_csv)
