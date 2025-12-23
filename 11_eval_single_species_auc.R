# Evaluate SINR-style model against IUCN range for a single species.
# Choose a target species by scientific name (must be in prototype_1000_species.csv)
# Load its IUCN range polygons from iucn_ranges_mammals_amphibians.rds (sf)
# Repair geometry and dissolve polygons (st_make_valid + summarise)
# Build a lon/lat grid over the IUCN bounding box (+ small margin)
# Sample random background points globally
# Extract WorldClim env variables with terra
# Build feature vectors (coord encoding + bioclims)
# Run SINR MLP to get probabilities for this species
# Use sf::st_intersects to label each point as inside/outside IUCN
# Compute AUC and Average Precision (AP)

library(dplyr)
library(readr)
library(torch)
library(terra)
library(pROC)
library(sf)

source("paths.R")
check_project_dirs()

# match "scientific_name" in prototype_1000_species.csv
target_species_name <- "Thylogale stigmatica"

# grid resolution (degrees) over IUCN bbox
grid_res <- 0.5

# number of random global background points
n_bg <- 2000L

# load metadata, model, species mapping

file_meta <- file.path(DIR_PROCESSED, "training_metadata_1000_species.rds")
if (!file.exists(file_meta)) {
  stop("Metadata file not found: ", file_meta,
       "\nRun 07_dataset_torch.R first.")
}
meta <- readRDS(file_meta)
feature_cols <- meta$feature_cols
input_dim    <- meta$input_dim
num_species  <- meta$num_species

message("Loaded metadata: input_dim = ", input_dim,
        ", num_species = ", num_species)

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
       "Check the spelling or choose a different species.")
}

species_label <- sp_row$prototype_index + 1L  # 1-based label index
message("Target species: ", target_species_name,
        " (label index = ", species_label, ")")

# load model architecture and weights
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

message("Loaded model weights from: ", model_path)

# load WorldClim env stack + IUCN polygons for this species

file_env_rds <- file.path(DIR_INTERMEDIATE, "worldclim21_bio_2-5m_stack.rds")
if (!file.exists(file_env_rds)) {
  stop("WorldClim env stack RDS not found: ", file_env_rds,
       "\nRun 02_prepare_worldclim.R first.")
}
env_stack <- readRDS(file_env_rds)

file_iucn_rds <- file.path(DIR_INTERMEDIATE, "iucn_ranges_mammals_amphibians.rds")
if (!file.exists(file_iucn_rds)) {
  stop("IUCN RDS file not found: ", file_iucn_rds,
       "\nCheck DIR_INTERMEDIATE or rerun 01_prepare_iucn.R.")
}

iucn_all <- readRDS(file_iucn_rds)

if (!inherits(iucn_all, "sf")) {
  stop("Expected iucn_ranges_mammals_amphibians.rds to be an sf object.\n",
       "Class is: ", paste(class(iucn_all), collapse = ", "))
}

if (!("sci_name" %in% names(iucn_all))) {
  stop("Could not find 'sci_name' column in IUCN sf object.\n",
       "Columns are: ", paste(names(iucn_all), collapse = ", "))
}

iucn_sp_sf <- iucn_all[iucn_all$sci_name == target_species_name, ]

if (nrow(iucn_sp_sf) == 0) {
  stop("No IUCN polygons found for species '", target_species_name, "'.")
}

message("Repairing IUCN polygon geometry for ", target_species_name, " ...")
iucn_sp_sf <- sf::st_make_valid(iucn_sp_sf)

#dissolve polygons into one MULTIPOLYGON to avoid slivers
iucn_sp_sf <- iucn_sp_sf %>% dplyr::summarise()

# Ensure CRS is WGS84
if (is.na(sf::st_crs(iucn_sp_sf))) {
  sf::st_crs(iucn_sp_sf) <- 4326
} else if (sf::st_crs(iucn_sp_sf)$epsg != 4326) {
  iucn_sp_sf <- sf::st_transform(iucn_sp_sf, 4326)
}

message("IUCN sf object ready with geometry type: ",
        unique(sf::st_geometry_type(iucn_sp_sf)))

# build a lon/lat grid over the IUCN bounding box

bbox <- sf::st_bbox(iucn_sp_sf)

margin <- 2  # degrees
lon_min <- max(as.numeric(bbox["xmin"]) - margin, -180)
lon_max <- min(as.numeric(bbox["xmax"]) + margin,  180)
lat_min <- max(as.numeric(bbox["ymin"]) - margin,  -90)
lat_max <- min(as.numeric(bbox["ymax"]) + margin,   90)

message("IUCN bounding box + margin:")
message("  lon: [", lon_min, ", ", lon_max, "]")
message("  lat: [", lat_min, ", ", lat_max, "]")

lon_seq <- seq(lon_min, lon_max, by = grid_res)
lat_seq <- seq(lat_min, lat_max, by = grid_res)

grid_df <- expand.grid(
  lon = lon_seq,
  lat = lat_seq
)

message("Grid (near range) has ", nrow(grid_df), " points.")

# Points for raster extraction (terra)
pts_grid_terra <- vect(
  data.frame(lon = grid_df$lon, lat = grid_df$lat),
  geom = c("lon", "lat"),
  crs  = "EPSG:4326"
)

# Points for polygon intersection (sf)
pts_grid_sf <- sf::st_as_sf(
  grid_df,
  coords = c("lon", "lat"),
  crs = 4326
)

# extract environment & features for grid

env_vals_grid <- terra::extract(env_stack, pts_grid_terra, df = TRUE)
env_vals_grid <- env_vals_grid[, -1, drop = FALSE]  # drop ID column
env_vals_grid[is.na(env_vals_grid)] <- 0

deg2rad <- function(x) x * pi / 180
lon_rad_grid <- deg2rad(sf::st_coordinates(pts_grid_sf)[,1])
lat_rad_grid <- deg2rad(sf::st_coordinates(pts_grid_sf)[,2])

feat_grid <- data.frame(
  x1 = sin(lon_rad_grid),
  x2 = cos(lon_rad_grid),
  x3 = sin(lat_rad_grid),
  x4 = cos(lat_rad_grid),
  env_vals_grid
)

available_feats <- intersect(feature_cols, names(feat_grid))
missing_feats   <- setdiff(feature_cols, names(feat_grid))

if (length(missing_feats) > 0) {
  warning("Some feature_cols are missing in feat_grid: ",
          paste(missing_feats, collapse = ", "),
          "\nUsing only the intersection of available features.")
}

feat_grid <- feat_grid[, available_feats, drop = FALSE]

complete_mask_grid <- complete.cases(feat_grid)
feat_grid <- feat_grid[complete_mask_grid, ]
pts_grid_sf <- pts_grid_sf[complete_mask_grid, ]

message("After cleaning, grid (near range) has ",
        nrow(feat_grid), " points.")

# label grid points as inside/outside IUCN using sf
inside_list_grid <- sf::st_intersects(pts_grid_sf, iucn_sp_sf, sparse = TRUE)
y_true_grid <- ifelse(lengths(inside_list_grid) > 0, 1L, 0L)

table_y_grid <- table(y_true_grid)
message("IUCN presence/absence on grid-only sample:")
print(table_y_grid)

# sample random global background points and build features

set.seed(123)
lon_bg <- runif(n_bg, min = -180, max = 180)
lat_bg <- runif(n_bg, min =  -90, max =  90)

bg_df <- data.frame(lon = lon_bg, lat = lat_bg)

# for raster extraction
pts_bg_terra <- vect(
  data.frame(lon = bg_df$lon, lat = bg_df$lat),
  geom = c("lon", "lat"),
  crs  = "EPSG:4326"
)

# for polygon intersection
pts_bg_sf <- sf::st_as_sf(
  bg_df,
  coords = c("lon", "lat"),
  crs = 4326
)

env_vals_bg <- terra::extract(env_stack, pts_bg_terra, df = TRUE)
env_vals_bg <- env_vals_bg[, -1, drop = FALSE]
env_vals_bg[is.na(env_vals_bg)] <- 0

lon_rad_bg <- deg2rad(sf::st_coordinates(pts_bg_sf)[,1])
lat_rad_bg <- deg2rad(sf::st_coordinates(pts_bg_sf)[,2])

feat_bg <- data.frame(
  x1 = sin(lon_rad_bg),
  x2 = cos(lon_rad_bg),
  x3 = sin(lat_rad_bg),
  x4 = cos(lat_rad_bg),
  env_vals_bg
)

feat_bg <- feat_bg[, available_feats, drop = FALSE]
complete_mask_bg <- complete.cases(feat_bg)
feat_bg <- feat_bg[complete_mask_bg, ]
pts_bg_sf <- pts_bg_sf[complete_mask_bg, ]

message("Background sample has ", nrow(feat_bg), " points after cleaning.")

# label background points as inside/outside IUCN using sf
inside_list_bg <- sf::st_intersects(pts_bg_sf, iucn_sp_sf, sparse = TRUE)
y_true_bg <- ifelse(lengths(inside_list_bg) > 0, 1L, 0L)

table_y_bg <- table(y_true_bg)
message("IUCN presence/absence on global background sample:")
print(table_y_bg)

# combine grid + background and run model

feat_all   <- rbind(feat_grid, feat_bg)
y_true_all <- c(y_true_grid, y_true_bg)

message("Combined sample size: ", length(y_true_all))

table_y_all <- table(y_true_all)
message("IUCN presence/absence on combined sample:")
print(table_y_all)

feat_mat_all <- as.matrix(feat_all[, available_feats, drop = FALSE])
x_all <- torch_tensor(feat_mat_all, dtype = torch_float(), device = device)

p_all <- model(x_all)
p_all_cpu <- p_all$to(device = "cpu")
p_all_mat <- as.matrix(as_array(p_all_cpu))

prob_species_all <- p_all_mat[, species_label]

message("Summary of predicted probabilities for ", target_species_name, ":")
print(summary(prob_species_all))

#compute AUC and average precision

if (length(unique(y_true_all)) < 2) {
  cat("\n[Warning] Combined y_true has only one class; ROC/AUC undefined.\n")
  auc_val <- NA_real_
} else {
  roc_obj <- roc(response = y_true_all, predictor = prob_species_all, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
}

ord <- order(prob_species_all, decreasing = TRUE)
y_sorted <- y_true_all[ord]

cum_tp    <- cumsum(y_sorted == 1L)
positions <- seq_along(y_sorted)

precisions <- cum_tp / positions

num_pos <- sum(y_sorted == 1L)
if (num_pos == 0) {
  ap_val <- NA_real_
} else {
  ap_val <- sum(precisions[y_sorted == 1L]) / num_pos
}

n_pos_all <- sum(y_true_all == 1L)
n_neg_all <- sum(y_true_all == 0L)

cat("\n===== Evaluation for", target_species_name, "=====\n")
cat("Number of points (combined):", length(y_true_all), "\n")
cat("Positives (inside IUCN):", n_pos_all, "\n")
cat("Negatives (outside IUCN):", n_neg_all, "\n")
cat(sprintf("AUC: %.4f\n", auc_val))
cat(sprintf("Average Precision (AP): %.4f\n", ap_val))
cat("=========================================\n")
