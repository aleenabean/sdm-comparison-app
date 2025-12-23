# use rgbif to download iNaturalist research-grade occurrences
# for the 1000 prototype species and save a cleaned occurrence table

library(dplyr)
library(readr)
library(rgbif)

source("paths.R")
check_project_dirs()


file_prototype <- file.path(DIR_INTERMEDIATE, "prototype_1000_species.csv")
if (!file.exists(file_prototype)) {
  stop("Prototype species list not found. Run 04_select_1000_species.R first.")
}

prototype_species <- read_csv(file_prototype, show_col_types = FALSE)

# Distinct GBIF taxonKeys for our 1000 species
taxon_keys <- prototype_species %>%
  distinct(gbif_usageKey) %>%
  pull(gbif_usageKey) %>%
  na.omit() %>%
  as.integer()

message("Number of distinct GBIF taxonKeys: ", length(taxon_keys))

if (length(taxon_keys) == 0) {
  stop("No GBIF taxonKeys found in prototype_1000_species.csv")
}

# submit GBIF occurrence download for iNaturalist dataset
# iNaturalist Research-grade Observations datasetKey on GBIF:
# https://www.gbif.org/dataset/50c9509d-22c7-4a22-a47d-8c48425ef4a7

inat_dataset_key <- "50c9509d-22c7-4a22-a47d-8c48425ef4a7"

message("Submitting GBIF occurrence download request...")

dw_key <- occ_download(
  pred("datasetKey", inat_dataset_key),
  pred_in("taxonKey", taxon_keys),
  pred("hasCoordinate", TRUE),
  format = "DWCA"
)

dw_key <- as.character(dw_key)
message("Download key: ", dw_key)
message("You can check status at: https://www.gbif.org/occurrence/download/", dw_key)

#get GBIF download
message("Waiting for GBIF download to finish (this can take a while)...")
occ_download_wait(dw_key)  

# occ_download_get() returns the downloaded zip path as a character string
zip_path <- occ_download_get(dw_key, overwrite = TRUE)
zip_path <- as.character(zip_path)

message("Download file path: ", zip_path)

# move the zip into rdata_raw/inat
if (!dir.exists(DIR_RAW_INAT)) {
  dir.create(DIR_RAW_INAT, recursive = TRUE)
}

dest_zip <- file.path(DIR_RAW_INAT, basename(zip_path))
file.rename(zip_path, dest_zip)
message("Moved GBIF DwC-A zip to: ", dest_zip)

#unzip and locate occurrence.txt

tmp_dir <- tempfile(pattern = "gbif_inat_dwca_")
dir.create(tmp_dir)
unzip(dest_zip, exdir = tmp_dir)

occ_file <- file.path(tmp_dir, "occurrence.txt")
if (!file.exists(occ_file)) {
  occ_candidates <- list.files(tmp_dir, pattern = "occurrence\\.txt$", full.names = TRUE, recursive = TRUE)
  if (length(occ_candidates) == 0) {
    stop("occurrence.txt not found in the DwC-A archive.")
  }
  occ_file <- occ_candidates[1]
}

message("Reading occurrence.txt from: ", occ_file)
message("This may take a while...")

occ <- read_tsv(
  occ_file,
  col_types = cols(.default = "c"),
  na = c("", "NA")
)

message("Raw occurrence rows: ", nrow(occ))

#clean and select columns we want
occ_clean <- occ %>%
  transmute(
    gbifID           = gbifID,
    taxonKey         = as.integer(taxonKey),
    scientificName   = scientificName,
    decimalLatitude  = as.numeric(decimalLatitude),
    decimalLongitude = as.numeric(decimalLongitude),
    eventDate        = eventDate
  ) %>%
  filter(
    !is.na(decimalLatitude),
    !is.na(decimalLongitude),
    !is.na(taxonKey)
  )

message("Cleaned occurrence rows (with coords & taxonKey): ", nrow(occ_clean))

# join to prototype species to attach prototype_index

prototype_minimal <- prototype_species %>%
  select(
    scientific_name,
    taxon_group,
    species_index,
    prototype_index,
    gbif_usageKey
  ) %>%
  rename(
    iucn_scientific_name = scientific_name,
    gbif_taxonKey        = gbif_usageKey
  )

occ_joined <- occ_clean %>%
  left_join(
    prototype_minimal,
    by = c("taxonKey" = "gbif_taxonKey")
  ) %>%
  filter(!is.na(prototype_index))

message("Occurrences matching prototype species: ", nrow(occ_joined))


# save cleaned occurrence table

file_occ <- file.path(DIR_INTERMEDIATE, "inat_occurrences_1000_species.csv")
write_csv(occ_joined, file_occ)

message("Saved cleaned iNat/GBIF occurrences to: ", file_occ)
