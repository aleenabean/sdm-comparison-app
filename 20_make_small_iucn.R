# create a smaller IUCN ranges file containing only the species used in the app.

library(dplyr)
library(readr)
library(sf)

source("paths.R")

#load the big IUCN RDS

big_iucn_path <- file.path(
  DIR_INTERMEDIATE,
  "iucn_ranges_mammals_amphibians_compressed.rds"
)

if (!file.exists(big_iucn_path)) {
  stop("Big IUCN RDS not found at: ", big_iucn_path)
}

cat("Reading big IUCN RDS (this may take a bit)...\n")
iucn_all <- readRDS(big_iucn_path)

if (!"sci_name" %in% names(iucn_all)) {
  stop("Expected column 'sci_name' in IUCN object but did not find it.\n",
       "Available columns: ", paste(names(iucn_all), collapse = ", "))
}

# get the 601 species actually used in the app and from SINR eval

eval_path <- file.path(DIR_PROCESSED, "sinr_1000species_iucn_eval.csv")
if (!file.exists(eval_path)) {
  stop("sinr_1000species_iucn_eval.csv not found at: ", eval_path)
}

sinr_eval <- readr::read_csv(eval_path, show_col_types = FALSE)

candidate_sp_cols <- c(
  "iucn_scientific_name",
  "scientific_name",
  "sci_name",
  "species",
  "species_name"
)

sp_hits <- intersect(candidate_sp_cols, names(sinr_eval))
if (length(sp_hits) == 0) {
  stop(
    "Could not find a species-name column in sinr_1000species_iucn_eval.csv.\n",
    "Looked for: ", paste(candidate_sp_cols, collapse = ", "), "\n",
    "Columns are: ", paste(names(sinr_eval), collapse = ", ")
  )
}

sp_col <- sp_hits[1]
cat("Using species column from eval file:", sp_col, "\n")

species_vec <- sinr_eval[[sp_col]]
species_choices <- sort(unique(species_vec))

cat("Number of app species (from SINR eval): ",
    length(species_choices), "\n")


cat("Number of species with non-NA sinr_label_index: ",
    length(species_choices), "\n")

# filter IUCN ranges down to those species

iucn_small <- iucn_all %>%
  filter(sci_name %in% species_choices)

cat("Original IUCN rows: ", nrow(iucn_all),
    " | Filtered rows: ", nrow(iucn_small), "\n")
#just check
if (nrow(iucn_small) == 0) {
  stop("After filtering, IUCN object has 0 rows – check that sci_name ",
       "matches iucn_scientific_name mapping.")
}

#save smaller, highly-compressed RDS

small_iucn_path <- file.path(
  DIR_INTERMEDIATE,
  "iucn_ranges_app_species.rds"
)

cat("Saving filtered IUCN object to:\n  ", small_iucn_path, "\n")
saveRDS(iucn_small, small_iucn_path, compress = "xz")

# show before/after sizes
big_size   <- file.info(big_iucn_path)$size / 1024^2
small_size <- file.info(small_iucn_path)$size / 1024^2

cat(sprintf("Big IUCN RDS size   : %.1f MB\n", big_size))
cat(sprintf("Small IUCN RDS size : %.1f MB\n", small_size))
cat("Done. Update the app to read iucn_ranges_app_species.rds instead.\n")
