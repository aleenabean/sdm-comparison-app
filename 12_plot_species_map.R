# plot SINR-predicted range map for a single species
# overlayed with its IUCN polygon

library(dplyr)
library(readr)
library(torch)
library(terra)
library(sf)

source("paths.R")
check_project_dirs()

target_species_name <- "Capra ibex"  # in prototype_1000_species.csv

# map grid resolution in degrees (finer than 0.5 for nicer maps)
grid_res <- 0.25

#load metadata, model, species mapping

file_meta <- file.path(DIR_PROCESSED, "training_metadata_1000_species.rds")
if (!file.exists(file_meta)) {
  stop("Metadata file not found: ", file_meta,
       "\nRun 07_dataset_torch.R first.")
}
meta <- readRDS(file_meta)
feature_cols <- meta$feature_cols
input_dim    <- meta$input_dim
num_species  <- meta$num_species

file_prototype <- file.path(DIR_INTERMEDIATE, "prototype_1000_species.csv")
if (!file.exists(file_prototype)) {
  stop("Prototype species list not found: ", file_prototype,
       "\nRun 04_select_1000_species.R first.")
}
prototype_species <- read_csv(file_prototype, show_col_types = FALSE)

sp_row <- prototype_species %>%
  filter(scientific_name == target_species_name) %>%
  slice(1)

if (nrow(sp_row) == 0) {
  stop("Target species '", target_species_name,
       "' not found in prototype_1000_species.csv.\n",
       "Check spelling or choose a different species.")
}

species_label <- sp_row$prototype_index + 1L

cat("Target species:", target_species_name,
    "(label index =", species_label, ")\n")

# load model
source("08_model_and_loss.R")

device <- torch_device("cpu")

model <- sinr_mlp(
  input_dim   = input_dim,
  num_species = num_species,
  hidden_dim  = 256L,
  dropout     = 0.2
)
model <- model$to(device = device)

model_path <- file.path(DIR_MODELS, "sinr_mlp_1000species_sinr.pt")
if (!file.exists(model_path)) {
  stop("Model weights not found: ", model_path,
       "\nRun 09_train_sinr_mlp.R first.")
}

state_dict <- torch_load(model_path)
model$load_state_dict(state_dict)
model$eval()

cat("Loaded model from:", model_path, "\n")

# load env stack + IUCN polygons (sf) and get bounding box

file_env_rds <- file.path(DIR_INTERMEDIATE, "worldclim21_bio_2-5m_stack.rds")
if (!file.exists(file_env_rds)) {
  stop("Env stack not found:", file_env_rds)
}
env_stack <- readRDS(file_env_rds)

file_iucn_rds <- file.path(DIR_INTERMEDIATE, "iucn_ranges_mammals_amphibians.rds")
if (!file.exists(file_iucn_rds)) {
  stop("IUCN RDS not found:", file_iucn_rds)
}

iucn_all <- readRDS(file_iucn_rds)

if (!inherits(iucn_all, "sf")) {
  stop("Expected IUCN RDS to be sf; got class:",
       paste(class(iucn_all), collapse = ", "))
}
if (!("sci_name" %in% names(iucn_all))) {
  stop("No 'sci_name' column in IUCN sf.\nCols: ",
       paste(names(iucn_all), collapse = ", "))
}

iucn_sp_sf <- iucn_all[iucn_all$sci_name == target_species_name, ]
if (nrow(iucn_sp_sf) == 0) {
  stop("No IUCN polygons found for species:", target_species_name)
}

# repair and dissolve geometry
iucn_sp_sf <- sf::st_make_valid(iucn_sp_sf)
iucn_sp_sf <- iucn_sp_sf %>% dplyr::summarise()

if (is.na(sf::st_crs(iucn_sp_sf))) {
  sf::st_crs(iucn_sp_sf) <- 4326
} else if (sf::st_crs(iucn_sp_sf)$epsg != 4326) {
  iucn_sp_sf <- sf::st_transform(iucn_sp_sf, 4326)
}

bbox <- sf::st_bbox(iucn_sp_sf)

# slight margin to show surroundings
margin <- 2
lon_min <- max(as.numeric(bbox["xmin"]) - margin, -180)
lon_max <- min(as.numeric(bbox["xmax"]) + margin,  180)
lat_min <- max(as.numeric(bbox["ymin"]) - margin,  -90)
lat_max <- min(as.numeric(bbox["ymax"]) + margin,   90)

cat("Map bounding box:\n")
cat("  lon:", lon_min, "to", lon_max, "\n")
cat("  lat:", lat_min, "to", lat_max, "\n")


# build grid over bbox + extract env + make features + run model


lon_seq <- seq(lon_min, lon_max, by = grid_res)
lat_seq <- seq(lat_min, lat_max, by = grid_res)

grid_df <- expand.grid(
  lon = lon_seq,
  lat = lat_seq
)

cat("Map grid has", nrow(grid_df), "cells.\n")

# for env extraction
pts_grid_terra <- vect(
  data.frame(lon = grid_df$lon, lat = grid_df$lat),
  geom = c("lon", "lat"),
  crs  = "EPSG:4326"
)

env_vals_grid <- terra::extract(env_stack, pts_grid_terra, df = TRUE)
env_vals_grid <- env_vals_grid[, -1, drop = TRUE]
env_vals_grid[is.na(env_vals_grid)] <- 0

# coord encoding
deg2rad <- function(x) x * pi / 180
lon_rad <- deg2rad(grid_df$lon)
lat_rad <- deg2rad(grid_df$lat)

feat_grid <- data.frame(
  x1 = sin(lon_rad),
  x2 = cos(lon_rad),
  x3 = sin(lat_rad),
  x4 = cos(lat_rad),
  env_vals_grid
)

available_feats <- intersect(feature_cols, names(feat_grid))
missing_feats   <- setdiff(feature_cols, names(feat_grid))
if (length(missing_feats) > 0) {
  warning("Missing features in map grid: ",
          paste(missing_feats, collapse = ", "),
          "\nUsing only available feature set.")
}
feat_grid <- feat_grid[, available_feats, drop = FALSE]

complete_mask <- complete.cases(feat_grid)
feat_grid <- feat_grid[complete_mask, ]
grid_df   <- grid_df[complete_mask, ]

cat("After cleaning, map grid has", nrow(grid_df), "cells.\n")

# run model
x_mat <- as.matrix(feat_grid)
x_t   <- torch_tensor(x_mat, dtype = torch_float(), device = device)

p_all <- model(x_t)
p_cpu <- p_all$to(device = "cpu")
p_mat <- as.matrix(as_array(p_cpu))

prob_species <- p_mat[, species_label]

cat("Predicted probability summary for", target_species_name, ":\n")
print(summary(prob_species))

# raster and plot with IUCN outline

# create SpatRaster from probabilities
rast_prob <- terra::rast(
  nrows = length(lat_seq),
  ncols = length(lon_seq),
  xmin  = lon_min,
  xmax  = lon_max,
  ymin  = lat_min,
  ymax  = lat_max,
  crs   = "EPSG:4326"
)

# Need to put probabilities into same ordering:
# expand.grid lat fastest / lon fastest can be tricky, so we reconstruct accordingly
# Our grid_df was made by expand.grid(lon, lat) => lon varies slowest, lat fastest.
# terra's rast is row-major in y (lat) then x (lon).
# We'll fill by matching each (lon,lat) to the correct cell index.

# Build a matrix [nrow(lat_seq), ncol(lon_seq)] in the right order
prob_mat <- matrix(NA_real_,
                   nrow = length(lat_seq),
                   ncol = length(lon_seq))

# map from lon/lat to indices
lon_index <- match(grid_df$lon, lon_seq)
lat_index <- match(grid_df$lat, lat_seq)

for (i in seq_along(prob_species)) {
  r <- lat_index[i]
  c <- lon_index[i]
  prob_mat[r, c] <- prob_species[i]
}

values(rast_prob) <- as.vector(prob_mat)

# convert IUCN sf to SpatVector for plotting
iucn_sp <- terra::vect(iucn_sp_sf)

# plot
plot(rast_prob,
     main = paste("SINR predicted probability for", target_species_name))
plot(iucn_sp, add = TRUE, border = "black", lwd = 2)
