# batch evaluation of SINR MLP against IUCN ranges for many species
# results are written to a CSV in DIR_EVAL

library(dplyr)
library(readr)
library(torch)
library(terra)
library(pROC)
library(sf)
library(purrr)

source("paths.R")
check_project_dirs()

if (!dir.exists(DIR_EVAL)) {
  dir.create(DIR_EVAL, recursive = TRUE)
}

num_species_to_eval <- 50L    # how many species to sample from prototype list
grid_res            <- 0.5    # degrees; coarser = faster
n_bg                <- 2000L  # number of global background points
min_pos_required    <- 5L     # skip AUC/AP if fewer positives than this
margin_deg          <- 2      # bbox margin in degrees

set.seed(123)  # for reproducible sampling / background points

output_csv <- file.path(DIR_EVAL, "sinr_1000species_iucn_eval.csv")

message("Will evaluate up to ", num_species_to_eval, " species.")
message("Results will be written to: ", output_csv)

# load metadata, model, env stack, IUCN ranges
#metadata
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

# prototype species list (names + label indices)
file_prototype <- file.path(DIR_INTERMEDIATE, "prototype_1000_species.csv")
if (!file.exists(file_prototype)) {
  stop("Prototype species list not found: ", file_prototype,
       "\nRun 04_select_1000_species.R first.")
}
prototype_species <- read_csv(file_prototype, show_col_types = FALSE)

# model
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

# env stack
file_env_rds <- file.path(DIR_INTERMEDIATE, "worldclim21_bio_2-5m_stack.rds")
if (!file.exists(file_env_rds)) {
  stop("Env stack not found: ", file_env_rds,
       "\nRun 02_prepare_worldclim.R first.")
}
env_stack <- readRDS(file_env_rds)

# IUCN ranges
file_iucn_rds <- file.path(DIR_INTERMEDIATE, "iucn_ranges_mammals_amphibians.rds")
if (!file.exists(file_iucn_rds)) {
  stop("IUCN RDS not found: ", file_iucn_rds,
       "\nRun 01_prepare_iucn.R first.")
}

iucn_all <- readRDS(file_iucn_rds)

if (!inherits(iucn_all, "sf")) {
  stop("Expected IUCN RDS to be sf; got class: ",
       paste(class(iucn_all), collapse = ", "))
}
if (!("sci_name" %in% names(iucn_all))) {
  stop("No 'sci_name' column in IUCN sf.\nCols: ",
       paste(names(iucn_all), collapse = ", "))
}


#choosing which species to evaluate
#only species present both in prototype list AND IUCN
proto_names <- prototype_species$scientific_name
iucn_names  <- unique(iucn_all$sci_name)

candidate_species <- intersect(proto_names, iucn_names)

if (length(candidate_species) == 0) {
  stop("No overlap between prototype species and IUCN sci_name.")
}

message("Number of species in both prototype & IUCN: ",
        length(candidate_species))

if (length(candidate_species) > num_species_to_eval) {
  eval_species <- sample(candidate_species, num_species_to_eval)
} else {
  eval_species <- candidate_species
}

message("Evaluating ", length(eval_species), " species.")

#get label index for a species
get_label_index <- function(sp_name) {
  row <- prototype_species %>%
    filter(scientific_name == sp_name) %>%
    slice(1)
  if (nrow(row) == 0) return(NA_integer_)
  row$prototype_index + 1L
}

#evaluate a single species

deg2rad <- function(x) x * pi / 180

eval_one_species <- function(sp_name) {
  message("\n--- Evaluating species: ", sp_name, " ---")
  
  label_idx <- get_label_index(sp_name)
  if (is.na(label_idx)) {
    warning("Species ", sp_name, " not in prototype_1000_species.csv; skipping.")
    return(tibble(
      scientific_name = sp_name,
      label_index     = NA_integer_,
      n_pos_grid      = NA_integer_,
      n_pos_bg        = NA_integer_,
      n_pos_all       = NA_integer_,
      n_neg_all       = NA_integer_,
      auc             = NA_real_,
      ap              = NA_real_,
      note            = "not_in_prototype"
    ))
  }
  
  # IUCN polygons for this species
  iucn_sp_sf <- iucn_all[iucn_all$sci_name == sp_name, ]
  if (nrow(iucn_sp_sf) == 0) {
    warning("No IUCN polygons found for ", sp_name, "; skipping.")
    return(tibble(
      scientific_name = sp_name,
      label_index     = label_idx,
      n_pos_grid      = NA_integer_,
      n_pos_bg        = NA_integer_,
      n_pos_all       = NA_integer_,
      n_neg_all       = NA_integer_,
      auc             = NA_real_,
      ap              = NA_real_,
      note            = "no_iucn_polygons"
    ))
  }
  
  # repair + dissolve geometry
  iucn_sp_sf <- sf::st_make_valid(iucn_sp_sf)
  iucn_sp_sf <- iucn_sp_sf %>% dplyr::summarise()
  
  if (is.na(sf::st_crs(iucn_sp_sf))) {
    sf::st_crs(iucn_sp_sf) <- 4326
  } else if (sf::st_crs(iucn_sp_sf)$epsg != 4326) {
    iucn_sp_sf <- sf::st_transform(iucn_sp_sf, 4326)
  }
  
  bbox <- sf::st_bbox(iucn_sp_sf)
  
  # build grid over bbox + margin
  lon_min <- max(as.numeric(bbox["xmin"]) - margin_deg, -180)
  lon_max <- min(as.numeric(bbox["xmax"]) + margin_deg,  180)
  lat_min <- max(as.numeric(bbox["ymin"]) - margin_deg,  -90)
  lat_max <- min(as.numeric(bbox["ymax"]) + margin_deg,   90)
  
  lon_seq <- seq(lon_min, lon_max, by = grid_res)
  lat_seq <- seq(lat_min, lat_max, by = grid_res)
  
  grid_df <- expand.grid(
    lon = lon_seq,
    lat = lat_seq
  )
  
  if (nrow(grid_df) == 0) {
    warning("Empty grid for ", sp_name, "; skipping.")
    return(tibble(
      scientific_name = sp_name,
      label_index     = label_idx,
      n_pos_grid      = NA_integer_,
      n_pos_bg        = NA_integer_,
      n_pos_all       = NA_integer_,
      n_neg_all       = NA_integer_,
      auc             = NA_real_,
      ap              = NA_real_,
      note            = "empty_grid"
    ))
  }
  
  # points for env extraction
  pts_grid_terra <- vect(
    data.frame(lon = grid_df$lon, lat = grid_df$lat),
    geom = c("lon", "lat"),
    crs  = "EPSG:4326"
  )
  
  # points for polygon intersection
  pts_grid_sf <- sf::st_as_sf(
    grid_df,
    coords = c("lon", "lat"),
    crs = 4326
  )
  
  # 3) extract env + build features for grid
  env_vals_grid <- terra::extract(env_stack, pts_grid_terra, df = TRUE)
  env_vals_grid <- env_vals_grid[, -1, drop = FALSE]
  env_vals_grid[is.na(env_vals_grid)] <- 0
  
  coords_grid <- sf::st_coordinates(pts_grid_sf)
  lon_rad_grid <- deg2rad(coords_grid[, 1])
  lat_rad_grid <- deg2rad(coords_grid[, 2])
  
  feat_grid <- data.frame(
    x1 = sin(lon_rad_grid),
    x2 = cos(lon_rad_grid),
    x3 = sin(lat_rad_grid),
    x4 = cos(lat_rad_grid),
    env_vals_grid
  )
  
  available_feats <- intersect(feature_cols, names(feat_grid))
  if (length(available_feats) == 0) {
    warning("No overlapping features for ", sp_name, "; skipping.")
    return(tibble(
      scientific_name = sp_name,
      label_index     = label_idx,
      n_pos_grid      = NA_integer_,
      n_pos_bg        = NA_integer_,
      n_pos_all       = NA_integer_,
      n_neg_all       = NA_integer_,
      auc             = NA_real_,
      ap              = NA_real_,
      note            = "no_matching_features"
    ))
  }
  
  feat_grid <- feat_grid[, available_feats, drop = FALSE]
  complete_mask_grid <- complete.cases(feat_grid)
  feat_grid <- feat_grid[complete_mask_grid, ]
  pts_grid_sf <- pts_grid_sf[complete_mask_grid, ]
  
  if (nrow(feat_grid) == 0) {
    warning("No complete grid features for ", sp_name, "; skipping.")
    return(tibble(
      scientific_name = sp_name,
      label_index     = label_idx,
      n_pos_grid      = NA_integer_,
      n_pos_bg        = NA_integer_,
      n_pos_all       = NA_integer_,
      n_neg_all       = NA_integer_,
      auc             = NA_real_,
      ap              = NA_real_,
      note            = "no_complete_grid_features"
    ))
  }
  
  # label grid points via sf::st_intersects
  inside_list_grid <- sf::st_intersects(pts_grid_sf, iucn_sp_sf, sparse = TRUE)
  y_true_grid <- ifelse(lengths(inside_list_grid) > 0, 1L, 0L)
  n_pos_grid <- sum(y_true_grid == 1L)
  
  #background
  lon_bg <- runif(n_bg, min = -180, max = 180)
  lat_bg <- runif(n_bg, min =  -90, max =  90)
  
  bg_df <- data.frame(lon = lon_bg, lat = lat_bg)
  
  pts_bg_terra <- vect(
    data.frame(lon = bg_df$lon, lat = bg_df$lat),
    geom = c("lon", "lat"),
    crs  = "EPSG:4326"
  )
  
  pts_bg_sf <- sf::st_as_sf(
    bg_df,
    coords = c("lon", "lat"),
    crs = 4326
  )
  
  env_vals_bg <- terra::extract(env_stack, pts_bg_terra, df = TRUE)
  env_vals_bg <- env_vals_bg[, -1, drop = FALSE]
  env_vals_bg[is.na(env_vals_bg)] <- 0
  
  coords_bg <- sf::st_coordinates(pts_bg_sf)
  lon_rad_bg <- deg2rad(coords_bg[, 1])
  lat_rad_bg <- deg2rad(coords_bg[, 2])
  
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
  
  inside_list_bg <- sf::st_intersects(pts_bg_sf, iucn_sp_sf, sparse = TRUE)
  y_true_bg <- ifelse(lengths(inside_list_bg) > 0, 1L, 0L)
  n_pos_bg <- sum(y_true_bg == 1L)
  
  # combine label vectors
  feat_all <- rbind(feat_grid, feat_bg)
  y_true_all <- c(y_true_grid, y_true_bg)
  
  n_pos_all <- sum(y_true_all == 1L)
  n_neg_all <- sum(y_true_all == 0L)
  
  note <- "ok"
  
  if (n_pos_all < min_pos_required) {
    warning("Species ", sp_name, " has only ", n_pos_all,
            " positives; skipping AUC/AP.")
    return(tibble(
      scientific_name = sp_name,
      label_index     = label_idx,
      n_pos_grid      = n_pos_grid,
      n_pos_bg        = n_pos_bg,
      n_pos_all       = n_pos_all,
      n_neg_all       = n_neg_all,
      auc             = NA_real_,
      ap              = NA_real_,
      note            = "too_few_positives"
    ))
  }
  
  if (length(unique(y_true_all)) < 2) {
    warning("Species ", sp_name, " has only one class in y_true_all.")
    return(tibble(
      scientific_name = sp_name,
      label_index     = label_idx,
      n_pos_grid      = n_pos_grid,
      n_pos_bg        = n_pos_bg,
      n_pos_all       = n_pos_all,
      n_neg_all       = n_neg_all,
      auc             = NA_real_,
      ap              = NA_real_,
      note            = "one_class_only"
    ))
  }
  
  # run model and compute AUC + AP
  x_mat <- as.matrix(feat_all)
  x_t   <- torch_tensor(x_mat, dtype = torch_float(), device = device)
  
  p_all <- model(x_t)
  p_all_cpu <- p_all$to(device = "cpu")
  p_all_mat <- as.matrix(as_array(p_all_cpu))
  
  prob_sp <- p_all_mat[, label_idx]
  
  #AUC
  roc_obj <- roc(response = y_true_all, predictor = prob_sp, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  
  #AP
  ord <- order(prob_sp, decreasing = TRUE)
  y_sorted <- y_true_all[ord]
  
  cum_tp    <- cumsum(y_sorted == 1L)
  positions <- seq_along(y_sorted)
  precisions <- cum_tp / positions
  
  num_pos <- sum(y_sorted == 1L)
  if (num_pos == 0) {
    ap_val <- NA_real_
    note <- paste0(note, ";no_positives_for_AP")
  } else {
    ap_val <- sum(precisions[y_sorted == 1L]) / num_pos
  }
  
  message("  n_pos_all = ", n_pos_all,
          ", n_neg_all = ", n_neg_all,
          ", AUC = ", sprintf("%.4f", auc_val),
          ", AP = ", sprintf("%.4f", ap_val))
  
  tibble(
    scientific_name = sp_name,
    label_index     = label_idx,
    n_pos_grid      = n_pos_grid,
    n_pos_bg        = n_pos_bg,
    n_pos_all       = n_pos_all,
    n_neg_all       = n_neg_all,
    auc             = auc_val,
    ap              = ap_val,
    note            = note
  )
}

#run evaluation loop

results_list <- vector("list", length(eval_species))

for (i in seq_along(eval_species)) {
  sp <- eval_species[i]
  message("\n===============================")
  message("Species ", i, " of ", length(eval_species), ": ", sp)
  message("===============================")
  
  res_i <- tryCatch(
    eval_one_species(sp),
    error = function(e) {
      warning("Error evaluating ", sp, ": ", conditionMessage(e))
      tibble(
        scientific_name = sp,
        label_index     = get_label_index(sp),
        n_pos_grid      = NA_integer_,
        n_pos_bg        = NA_integer_,
        n_pos_all       = NA_integer_,
        n_neg_all       = NA_integer_,
        auc             = NA_real_,
        ap              = NA_real_,
        note            = paste0("error: ", conditionMessage(e))
      )
    }
  )
  
  results_list[[i]] <- res_i
}

results <- bind_rows(results_list)

message("\nFinished batch evaluation for ", nrow(results), " species.")
message("Writing results to: ", output_csv)

write_csv(results, output_csv)

message("Done.")
