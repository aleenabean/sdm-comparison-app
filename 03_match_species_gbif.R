# Use rgbif to match a subset of IUCN species names to GBIF taxon keys.

library(dplyr)
library(readr)
library(purrr)
library(rgbif)

source("paths.R")
check_project_dirs()

file_species_csv <- file.path(DIR_INTERMEDIATE, "iucn_species_table_mammals_amphibians.csv")
if (!file.exists(file_species_csv)) {
  stop("IUCN species table not found. Run 01_prepare_iucn.R first.")
}

iucn_species <- read_csv(file_species_csv, show_col_types = FALSE)

#choosing a subset for a 1000-species prototype
#sample 15000 species
set.seed(42)
target_species <- iucn_species %>%
  arrange(taxon_group, scientific_name) %>%
  distinct(scientific_name, taxon_group, species_index) %>%
  slice_sample(n = 1500)

message("Number of species in subset: ", nrow(target_species))


lookup_gbif <- function(name) {
  res <- tryCatch(
    {
      rgbif::name_backbone(name = name)
    },
    error = function(e) {
      message("GBIF lookup failed for: ", name, " (", conditionMessage(e), ")")
      return(NULL)
    }
  )
  
  if (is.null(res)) {
    tibble(
      scientific_name = name,
      gbif_usageKey = NA_integer_,
      gbif_rank = NA_character_,
      gbif_status = NA_character_
    )
  } else {
    tibble(
      scientific_name = name,
      gbif_usageKey = res$usageKey,
      gbif_rank = res$rank,
      gbif_status = res$status
    )
  }
}


gbif_matches <- target_species %>%
  distinct(scientific_name) %>%
  mutate(row_id = row_number()) %>%
  group_by(row_id) %>%
  reframe(
    lookup_gbif(scientific_name)
  ) %>%
  select(-row_id)

species_gbif_map_subset <- target_species %>%
  left_join(gbif_matches, by = "scientific_name")

file_species_gbif <- file.path(DIR_INTERMEDIATE, "iucn_species_gbif_map_1500.csv")
write_csv(species_gbif_map_subset, file_species_gbif)

message("Saved subset IUCN↔GBIF species map to: ", file_species_gbif)

#summary of how many matched
summary_tbl <- species_gbif_map_subset %>%
  distinct(scientific_name, gbif_usageKey, gbif_status) %>%
  mutate(has_match = !is.na(gbif_usageKey)) %>%
  count(has_match)

print(summary_tbl)
