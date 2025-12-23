# Build a training-ready table:
# *cap observations per species
# *attach WorldClim env covariates
# *add coord encoding
# *create train/val/test split

library(dplyr)
library(readr)
library(terra)

source("paths.R")
check_project_dirs()


# get occurrences and WorldClim 

file_occ <- file.path(DIR_INTERMEDIATE, "inat_occurrences_1000_species.csv")
if (!file.exists(file_occ)) {
  stop("Occurrence file not found: ", file_occ, 
       "\nRun 05_download_inat_occurrences.R first.")
}

occ <- read_csv(file_occ, show_col_types = FALSE)
message("Loaded occurrences: ", nrow(occ))

file_env_rds <- file.path(DIR_INTERMEDIATE, "worldclim21_bio_2-5m_stack.rds")
if (!file.exists(file_env_rds)) {
  stop("WorldClim env stack RDS not found: ", file_env_rds,
       "\nRun 02_prepare_worldclim.R first.")
}

env_stack <- readRDS(file_env_rds)
message("Loaded WorldClim stack with ", nlyr(env_stack), " layers.")


# cap observations per species at 1000 (prototype_index)
cap_n <- 1000L

occ_capped <- occ %>%
  filter(!is.na(prototype_index)) %>%
  group_by(prototype_index) %>%
  group_modify(
    ~ {
      n_rows <- nrow(.x)
      slice_sample(.x, n = min(cap_n, n_rows))
    }
  ) %>%
  ungroup()

message("Occurrences after capping at ", cap_n, " per species: ", 
        nrow(occ_capped))

# build SpatVector of points for extraction

pts <- vect(
  occ_capped %>%
    transmute(
      lon = decimalLongitude,
      lat = decimalLatitude
    ),
  geom = c("lon", "lat"),
  crs  = "EPSG:4326"
)

message("Created SpatVector with ", nrow(occ_capped), " points.")

# extract WorldClim covariates for each point

message("Extracting WorldClim covariates at point locations (this may take a bit)...")

env_vals <- terra::extract(env_stack, pts, df = TRUE)


env_vals <- env_vals[, -1, drop = FALSE]

message("Extracted env matrix with ", nrow(env_vals), " rows and ", 
        ncol(env_vals), " env variables.")

if (nrow(env_vals) != nrow(occ_capped)) {
  stop("Row mismatch between occurrences and extracted env values.")
}

# combine occurrences + env covariates

env_names <- names(env_stack)
names(env_vals) <- env_names

occ_env <- bind_cols(
  occ_capped %>%
    mutate(row_id = row_number()),
  env_vals
)

message("Combined occurrences + env vars -> rows: ", nrow(occ_env), 
        ", columns: ", ncol(occ_env))

# adding 4D coordinate encoding

deg2rad <- function(x) x * pi / 180

occ_env <- occ_env %>%
  mutate(
    lon_rad = deg2rad(decimalLongitude),
    lat_rad = deg2rad(decimalLatitude),
    x1 = sin(lon_rad),
    x2 = cos(lon_rad),
    x3 = sin(lat_rad),
    x4 = cos(lat_rad)
  )

#train/val/test split per species

set.seed(42)

occ_split <- occ_env %>%
  group_by(prototype_index) %>%
  mutate(
    .u = runif(dplyr::n()),
    split = case_when(
      .u < 0.8 ~ "train",
      .u < 0.9 ~ "val",
      TRUE     ~ "test"
    )
  ) %>%
  ungroup() %>%
  select(-.u)

table_splits <- occ_split %>%
  count(split)

message("Split counts:")
print(table_splits)

# save processed table

if (!dir.exists(DIR_PROCESSED)) {
  dir.create(DIR_PROCESSED, recursive = TRUE)
}

file_processed <- file.path(DIR_PROCESSED, "training_points_1000_species.csv")
write_csv(occ_split, file_processed)

message("Saved processed training table to: ", file_processed)
