# get 1000 prototype species from the 1500 IUCN/GBIF subset

library(dplyr)
library(readr)

source("paths.R")
check_project_dirs()

file_species_gbif <- file.path(DIR_INTERMEDIATE, "iucn_species_gbif_map_1500.csv")
if (!file.exists(file_species_gbif)) {
  stop("Subset species map not found. Run 03_match_species_gbif.R first.")
}

species_gbif <- read_csv(file_species_gbif, show_col_types = FALSE)

# keeping only species with a GBIF match
matched_species <- species_gbif %>%
  filter(!is.na(gbif_usageKey)) %>%
  distinct(scientific_name, taxon_group, species_index, gbif_usageKey, gbif_rank, gbif_status)

message("Species with GBIF match: ", nrow(matched_species))

#want to balance mammals and amphibians
set.seed(42)
matched_species %>% count(taxon_group) %>% print()
total_target <- 1000
group_counts <- matched_species %>%
  count(taxon_group, name = "n_group") %>%
  mutate(
    frac       = n_group / sum(n_group),
    target_n   = round(frac * total_target)
  )

print(group_counts)

selected_list <- group_counts %>%
  rowwise() %>%
  do({
    tg  <- .$taxon_group
    n_t <- .$target_n
    
    pool <- matched_species %>% filter(taxon_group == tg)
    n_take <- min(n_t, nrow(pool))
    
    pool %>%
      slice_sample(n = n_take)
  })

prototype_species <- selected_list %>%
  ungroup() %>%
  arrange(taxon_group, scientific_name) %>%
  mutate(
    
    prototype_index = row_number() - 1L
  )

file_prototype <- file.path(DIR_INTERMEDIATE, "prototype_1000_species.csv")
write_csv(prototype_species, file_prototype)

message("Saved prototype 1000-species list to: ", file_prototype)
