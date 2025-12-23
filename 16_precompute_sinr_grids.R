# ran LOCALLY (not on shinyapps.io)
# loads the trained SINR model (using torch)
# loops over all species used in the app (species_choices)
# predicts SINR on a grid around the IUCN range
# saves one small RDS per species in data_processed/

# Shiny app can use these RDS files to remove torch 

library(dplyr)
library(readr)
library(terra)
library(sf)
library(torch)

source("paths.R")
check_project_dirs()


# worldClim env stack (already saved as a SpatRaster RDS)
FILE_ENV_RDS  <- file.path(DIR_INTERMEDIATE, "worldclim21_bio_2-5m_stack.rds")
if (!file.exists(FILE_ENV_RDS)) {
  stop("Env stack RDS not found at: ", FILE_ENV_RDS,
       "\nRun 02_prepare_worldclim.R first.")
}
env_stack <- readRDS(FILE_ENV_RDS)

# iNat / GBIF occurrences
FILE_INAT_CSV <- file.path(DIR_INTERMEDIATE, "inat_occurrences_1000_species.csv")
if (!file.exists(FILE_INAT_CSV)) {
  stop("inat_occurrences_1000_species.csv not found at: ", FILE_INAT_CSV)
}
inat_occ <- readr::read_csv(FILE_INAT_CSV, show_col_types = FALSE)

# IUCN ranges (BIG file, must exist locally for this script)
FILE_IUCN_RDS <- file.path(DIR_INTERMEDIATE, "iucn_ranges_mammals_amphibians_compressed.rds")
if (!file.exists(FILE_IUCN_RDS)) {
  stop("IUCN ranges RDS not found at: ", FILE_IUCN_RDS,
       "\nPut iucn_ranges_mammals_amphibians_compressed.rds back into data_intermediate/.")
}
iucn_all <- readRDS(FILE_IUCN_RDS)

# prototype species + training metadata
FILE_PROTO <- file.path(DIR_INTERMEDIATE, "prototype_1000_species.csv")
if (!file.exists(FILE_PROTO)) {
  stop("prototype_1000_species.csv not found at: ", FILE_PROTO)
}
prototype_species <- readr::read_csv(FILE_PROTO, show_col_types = FALSE)

FILE_META <- file.path(DIR_PROCESSED, "training_metadata_1000_species.rds")
if (!file.exists(FILE_META)) {
  stop("training_metadata_1000_species.rds not found at: ", FILE_META)
}
meta         <- readRDS(FILE_META)
feature_cols <- meta$feature_cols
input_dim    <- meta$input_dim
num_species  <- meta$num_species


# iNat species column
sp_col_inat <- "iucn_scientific_name"
if (!sp_col_inat %in% names(inat_occ)) {
  stop("Column 'iucn_scientific_name' not found in inat_occ.")
}

#species in each source
inat_spp  <- unique(na.omit(inat_occ[[sp_col_inat]]))
proto_spp <- prototype_species$scientific_name
iucn_spp  <- unique(iucn_all$sci_name)

# final list used in app
species_choices <- sort(intersect(intersect(inat_spp, proto_spp), iucn_spp))

length(species_choices)
head(species_choices, 10)


#helping stuff

GRID_RES_GLOBAL <- 0.25
deg2rad <- function(x) x * pi / 180

get_label_index <- function(sp_name) {
  row <- prototype_species %>%
    dplyr::filter(scientific_name == sp_name) %>%
    dplyr::slice(1)
  
  if (nrow(row) == 0) {
    return(NA_integer_)
  }
  
  # torch label is prototype_index + 1
  as.integer(row$prototype_index[1]) + 1L
}


get_iucn_polygon <- function(sp_name) {
  iucn_sp_sf <- iucn_all[iucn_all$sci_name == sp_name, ]
  if (nrow(iucn_sp_sf) == 0) return(NULL)
  iucn_sp_sf <- sf::st_make_valid(iucn_sp_sf)
  iucn_sp_sf <- iucn_sp_sf %>% dplyr::summarise()
  if (is.na(sf::st_crs(iucn_sp_sf))) {
    sf::st_crs(iucn_sp_sf) <- 4326
  } else if (sf::st_crs(iucn_sp_sf)$epsg != 4326) {
    iucn_sp_sf <- sf::st_transform(iucn_sp_sf, 4326)
  }
  iucn_sp_sf
}

get_bbox_for_species <- function(sp_name, margin = 1) {
  iucn_sp_sf <- get_iucn_polygon(sp_name)
  if (!is.null(iucn_sp_sf)) {
    bb <- sf::st_bbox(iucn_sp_sf)
    c(
      max(as.numeric(bb["xmin"]) - margin, -180),
      min(as.numeric(bb["xmax"]) + margin,  180),
      max(as.numeric(bb["ymin"]) - margin,  -90),
      min(as.numeric(bb["ymax"]) + margin,   90)
    )
  } else {
    c(-20, 20, -20, 20)
  }
}

slugify_species <- function(sp_name) {
  s <- tolower(sp_name)
  s <- gsub("[^a-z0-9]+", "_", s)
  s <- gsub("^_|_$", "", s)
  s
}


sinr_mlp <- nn_module(
  "sinr_mlp",
  initialize = function(input_dim,
                        num_species,
                        hidden_dim = 256L,
                        dropout = 0.2) {
    self$fc1 <- nn_linear(input_dim, hidden_dim)
    self$fc2 <- nn_linear(hidden_dim, hidden_dim)
    self$fc3 <- nn_linear(hidden_dim, hidden_dim)
    self$out <- nn_linear(hidden_dim, num_species)
    self$dropout_p <- dropout
  },
  forward = function(x) {
    h1 <- self$fc1(x)
    h1 <- nnf_relu(h1)
    h2 <- self$fc2(h1)
    h2 <- nnf_relu(h2)
    h2 <- h2 + h1
    h3 <- self$fc3(h2)
    h3 <- nnf_relu(h3)
    h3 <- nnf_dropout(h3, p = self$dropout_p, training = self$training)
    logits <- self$out(h3)
    logits$sigmoid()
  }
)

device <- torch_device("cpu")

sinr_model <- sinr_mlp(
  input_dim   = input_dim,
  num_species = num_species,
  hidden_dim  = 256L,
  dropout     = 0.2
)
sinr_model$to(device = device)

MODEL_PATH <- file.path(DIR_MODELS, "sinr_mlp_1000species_sinr.pt")
if (!file.exists(MODEL_PATH)) {
  stop("SINR model file not found at: ", MODEL_PATH)
}
state_dict <- torch_load(MODEL_PATH)
sinr_model$load_state_dict(state_dict)
sinr_model$eval()

#single-species SINR grid prediction

predict_sinr_to_grid_once <- function(sp_name, bbox) {
  label_idx <- get_label_index(sp_name)
  if (is.na(label_idx)) {
    stop("Species not in prototype_1000_species.csv: ", sp_name)
  }
  
  lon_min <- bbox[1]; lon_max <- bbox[2]
  lat_min <- bbox[3]; lat_max <- bbox[4]
  lon_min <- max(lon_min, -180); lon_max <- min(lon_max, 180)
  lat_min <- max(lat_min,  -90); lat_max <- min(lat_max, 90)
  
  lon_seq <- seq(lon_min, lon_max, by = GRID_RES_GLOBAL)
  lat_seq <- seq(lat_min, lat_max, by = GRID_RES_GLOBAL)
  grid_df <- expand.grid(lon = lon_seq, lat = lat_seq)
  if (nrow(grid_df) == 0) stop("Empty grid for ", sp_name)
  
  pts_grid <- terra::vect(grid_df, geom = c("lon", "lat"), crs = "EPSG:4326")
  env_vals <- terra::extract(env_stack, pts_grid, df = TRUE)
  env_vals <- env_vals[, -1, drop = FALSE]
  env_vals[is.na(env_vals)] <- 0
  
  lon_rad <- deg2rad(grid_df$lon)
  lat_rad <- deg2rad(grid_df$lat)
  
  feat_grid <- data.frame(
    x1 = sin(lon_rad),
    x2 = cos(lon_rad),
    x3 = sin(lat_rad),
    x4 = cos(lat_rad),
    env_vals
  )
  
  available_feats <- intersect(feature_cols, names(feat_grid))
  feat_grid <- feat_grid[, available_feats, drop = FALSE]
  
  complete_mask <- stats::complete.cases(feat_grid)
  feat_grid <- feat_grid[complete_mask, ]
  grid_df   <- grid_df[complete_mask, ]
  
  # chunked prediction to avoid huge as.matrix()
  n_rows <- nrow(feat_grid)
  prob_sp <- numeric(n_rows)
  
  chunk_size <- 50000L  
  
  starts <- seq(1L, n_rows, by = chunk_size)
  for (s in starts) {
    e <- min(s + chunk_size - 1L, n_rows)
    idx <- s:e
    
    mat_chunk <- as.matrix(feat_grid[idx, , drop = FALSE])
    
    x_t   <- torch_tensor(mat_chunk, dtype = torch_float(), device = device)
    p_all <- sinr_model(x_t)
    p_cpu <- p_all$to(device = "cpu")
    p_mat <- as.matrix(as_array(p_cpu))
    
    prob_sp[idx] <- p_mat[, label_idx]
    
    rm(mat_chunk, x_t, p_all, p_cpu, p_mat)
    gc(verbose = FALSE)
  }
  
  grid_df$prob <- prob_sp
  
  n_lon <- length(lon_seq)
  n_lat <- length(lat_seq)
  prob_mat <- matrix(NA_real_, nrow = n_lat, ncol = n_lon)
  lon_index <- match(grid_df$lon, lon_seq)
  lat_index <- match(grid_df$lat, lat_seq)
  for (i in seq_along(prob_sp)) {
    r <- lat_index[i]; c <- lon_index[i]
    prob_mat[r, c] <- prob_sp[i]
  }
  
  rast_prob <- terra::rast(
    nrows = n_lat,
    ncols = n_lon,
    xmin  = lon_min,
    xmax  = lon_max,
    ymin  = lat_min,
    ymax  = lat_max,
    vals  = as.vector(t(prob_mat)),
    crs   = "EPSG:4326"
  )
  
  list(
    raster  = rast_prob,
    grid_df = grid_df
  )
}

#loop over all species and save RDS files

dir.create(DIR_PROCESSED, showWarnings = FALSE, recursive = TRUE)

for (sp in species_choices) {
  cat("Precomputing SINR grid for:", sp, "...\n")
  bbox <- get_bbox_for_species(sp, margin = 1)
  proj <- predict_sinr_to_grid_once(sp, bbox)
  fname <- file.path(DIR_PROCESSED, paste0("sinr_grid_", slugify_species(sp), ".rds"))
  saveRDS(proj, fname)
}

cat("Done. Precomputed SINR grids written to: ", DIR_PROCESSED, "\n")
