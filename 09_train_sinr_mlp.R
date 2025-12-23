# train SINR-style MLP with single-positive multi-label loss
# and random background locations

library(dplyr)
library(torch)
library(terra)
library(coro)
library(readr)

source("paths.R")
check_project_dirs()

#load datasets / dataloaders / metadata + model & loss

source("07_dataset_torch.R")   # defines train_dl, val_dl, train_data, feature_cols, input_dim, num_species
source("08_model_and_loss.R")  # defines sinr_mlp, sinr_spl_loss

#load env stack for random features
file_env_rds <- file.path(DIR_INTERMEDIATE, "worldclim21_bio_2-5m_stack.rds")
if (!file.exists(file_env_rds)) {
  stop("WorldClim env stack RDS not found: ", file_env_rds,
       "\nRun 02_prepare_worldclim.R first.")
}
env_stack <- readRDS(file_env_rds)

env_names <- names(env_stack)
message("Training with input_dim = ", input_dim, 
        ", num_species = ", num_species)


device <- torch_device("cpu")

set.seed(42)
torch_manual_seed(42)

#sample random-location features

sample_random_features <- function(batch_size, env_stack, feature_cols) {
  lon <- runif(batch_size, min = -180, max = 180)
  lat <- runif(batch_size, min = -90,  max = 90)
  
  pts <- vect(
    data.frame(lon = lon, lat = lat),
    geom = c("lon", "lat"),
    crs  = "EPSG:4326"
  )
  
  env_vals <- terra::extract(env_stack, pts, df = TRUE)
  env_vals <- env_vals[, -1, drop = FALSE]  # drop ID
  
  env_vals[is.na(env_vals)] <- 0
  
  deg2rad <- function(x) x * pi / 180
  lon_rad <- deg2rad(lon)
  lat_rad <- deg2rad(lat)
  
  df_rand <- data.frame(
    x1 = sin(lon_rad),
    x2 = cos(lon_rad),
    x3 = sin(lat_rad),
    x4 = cos(lat_rad),
    env_vals
  )
  
  # keeping only features we actually use
  available_feats <- intersect(feature_cols, names(df_rand))
  if (length(available_feats) != length(feature_cols)) {
    warning("Some feature_cols missing in random features; using intersection only.")
  }
  
  as.matrix(df_rand[, available_feats, drop = FALSE])
}

#model & optimizer

model <- sinr_mlp(
  input_dim   = input_dim,
  num_species = num_species,
  hidden_dim  = 256L,
  dropout     = 0.2
)

model <- model$to(device = device)

learning_rate <- 1e-3
optimizer <- optim_adam(model$parameters, lr = learning_rate)

lambda_pos <- 50    # positive weight in loss
num_epochs <- 10L 

#training & validation loops

train_one_epoch <- function(epoch,
                            model,
                            train_dl,
                            optimizer,
                            env_stack,
                            feature_cols,
                            lambda_pos,
                            device) {
  model$train()
  
  total_loss <- 0
  batch_count <- 0
  
  coro::loop(for (batch in train_dl) {
    optimizer$zero_grad()
    
    x <- batch$x$to(device = device)
    y <- batch$y$to(device = device)
    
    p_data <- model(x)  
    
    # random background batch (same size as data)
    B <- x$size(1)
    rand_mat <- sample_random_features(
      batch_size   = as.integer(B),
      env_stack    = env_stack,
      feature_cols = feature_cols
    )
    x_rand <- torch_tensor(rand_mat, dtype = torch_float(), device = device)
    
    p_rand <- model(x_rand) 
    
    # SINR-style loss
    loss <- sinr_spl_loss(
      p_data   = p_data,
      y        = y,
      p_rand   = p_rand,
      lambda_pos = lambda_pos
    )
    
    loss$backward()
    optimizer$step()
    
    total_loss <- total_loss + as.numeric(loss$item())
    batch_count <- batch_count + 1L
    
    if (batch_count %% 50 == 0) {
      cat(sprintf("  [Epoch %d] Batch %d, loss = %.4f\n",
                  epoch, batch_count, as.numeric(loss$item())))
    }
  })
  
  mean_loss <- total_loss / max(1, batch_count)
  cat(sprintf("Epoch %d: mean train loss = %.4f\n", epoch, mean_loss))
  
  mean_loss
}

evaluate_on_val <- function(model,
                            val_dl,
                            env_stack,
                            feature_cols,
                            lambda_pos,
                            device) {
  model$eval()
  
  total_loss <- 0
  batch_count <- 0
  
  coro::loop(for (batch in val_dl) {
    x <- batch$x$to(device = device)
    y <- batch$y$to(device = device)
    
    p_data <- model(x)
    
    B <- x$size(1)
    rand_mat <- sample_random_features(
      batch_size   = as.integer(B),
      env_stack    = env_stack,
      feature_cols = feature_cols
    )
    x_rand <- torch_tensor(rand_mat, dtype = torch_float(), device = device)
    p_rand <- model(x_rand)
    
    loss <- sinr_spl_loss(
      p_data   = p_data,
      y        = y,
      p_rand   = p_rand,
      lambda_pos = lambda_pos
    )
    
    total_loss <- total_loss + as.numeric(loss$item())
    batch_count <- batch_count + 1L
  })
  
  mean_loss <- total_loss / max(1, batch_count)
  cat(sprintf("Validation mean loss = %.4f\n", mean_loss))
  
  mean_loss
}

#running training

train_history <- data.frame(
  epoch = integer(),
  train_loss = numeric(),
  val_loss = numeric()
)

for (epoch in seq_len(num_epochs)) {
  cat("=== Epoch", epoch, "===\n")
  
  train_loss <- train_one_epoch(
    epoch       = epoch,
    model       = model,
    train_dl    = train_dl,
    optimizer   = optimizer,
    env_stack   = env_stack,
    feature_cols = feature_cols,
    lambda_pos  = lambda_pos,
    device      = device
  )
  
  val_loss <- evaluate_on_val(
    model       = model,
    val_dl      = val_dl,
    env_stack   = env_stack,
    feature_cols = feature_cols,
    lambda_pos  = lambda_pos,
    device      = device
  )
  
  train_history <- rbind(
    train_history,
    data.frame(
      epoch = epoch,
      train_loss = train_loss,
      val_loss = val_loss
    )
  )
}

#save model & history

if (!dir.exists(DIR_MODELS)) {
  dir.create(DIR_MODELS, recursive = TRUE)
}
if (!dir.exists(DIR_EVAL)) {
  dir.create(DIR_EVAL, recursive = TRUE)
}

model_path <- file.path(DIR_MODELS, "sinr_mlp_1000species_sinr.pt")
torch_save(model$state_dict(), model_path)
cat("Saved SINR-style model weights to: ", model_path, "\n")

history_path <- file.path(DIR_EVAL, "train_history_1000species_sinr.csv")
write_csv(train_history, history_path)
cat("Saved training history to: ", history_path, "\n")
