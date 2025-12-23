# build R+torch datasets and dataloaders for the 1000-species prototype

library(dplyr)
library(readr)
library(torch)
library(terra)

source("paths.R")
check_project_dirs()

# load processed training table and WorldClim stack

file_processed <- file.path(DIR_PROCESSED, "training_points_1000_species.csv")
if (!file.exists(file_processed)) {
  stop("Processed training file not found: ", file_processed,
       "\nRun 06_build_training_table.R first.")
}

data_all <- read_csv(file_processed, show_col_types = FALSE)
message("Loaded processed training data: ", nrow(data_all), " rows.")

file_env_rds <- file.path(DIR_INTERMEDIATE, "worldclim21_bio_2-5m_stack.rds")
if (!file.exists(file_env_rds)) {
  stop("WorldClim env stack RDS not found: ", file_env_rds,
       "\nRun 02_prepare_worldclim.R first.")
}
env_stack <- readRDS(file_env_rds)
env_names <- names(env_stack)

message("Env variable names from WorldClim:")
print(env_names)

#define feature + label columns

coord_cols <- c("x1", "x2", "x3", "x4")

missing_coord <- setdiff(coord_cols, names(data_all))
if (length(missing_coord) > 0) {
  stop("Missing coord encoding columns in data: ", paste(missing_coord, collapse = ", "))
}

missing_env <- setdiff(env_names, names(data_all))
if (length(missing_env) > 0) {
  warning(
    "Some env variables from env_stack are missing in training data:\n",
    paste(missing_env, collapse = ", "),
    "\nWe will only use those that are present."
  )
  env_used <- intersect(env_names, names(data_all))
} else {
  env_used <- env_names
}

feature_cols <- c(coord_cols, env_used)

message("Number of feature columns: ", length(feature_cols))
message("Feature columns:")
print(feature_cols)

if (!("prototype_index" %in% names(data_all))) {
  stop("prototype_index column not found in training data.")
}

# 1-based labels for torch
data_all <- data_all %>%
  mutate(label = prototype_index + 1L)

#get rid of incomplete feature rows

complete_mask <- stats::complete.cases(data_all[, feature_cols])
data_all_clean <- data_all[complete_mask, ]

message("Rows after dropping incomplete feature cases: ", nrow(data_all_clean))

# train / val / test splits

if (!("split" %in% names(data_all_clean))) {
  stop("split column not found in training data (expected 'train', 'val', 'test').")
}

split_counts <- data_all_clean %>%
  count(split)

message("Split counts after cleaning:")
print(split_counts)

train_data <- data_all_clean %>% filter(split == "train")
val_data   <- data_all_clean %>% filter(split == "val")
test_data  <- data_all_clean %>% filter(split == "test")

message("Train rows: ", nrow(train_data))
message("Val rows:   ", nrow(val_data))
message("Test rows:  ", nrow(test_data))

# torch dataset class

species_dataset <- dataset(
  name = "SpeciesDataset",
  
  initialize = function(df, feature_cols, label_col) {
    self$df <- df
    self$feature_cols <- feature_cols
    self$label_col <- label_col
    
    self$num_features <- length(feature_cols)
    self$num_rows <- nrow(df)
  },
  
  .getitem = function(i) {
    row <- self$df[i, , drop = FALSE]
    x_vec <- as.numeric(row[1, self$feature_cols])
    
    x <- torch_tensor(x_vec, dtype = torch_float())
    y <- torch_tensor(as.integer(row[[self$label_col]]), dtype = torch_long())
    
    list(x = x, y = y)
  },
  
  .length = function() {
    self$num_rows
  }
)

#datasets and dataloaders

train_dataset <- species_dataset(
  df = train_data,
  feature_cols = feature_cols,
  label_col = "label"
)

val_dataset <- species_dataset(
  df = val_data,
  feature_cols = feature_cols,
  label_col = "label"
)

batch_size_train <- 512L
batch_size_val   <- 1024L

train_dl <- dataloader(
  train_dataset,
  batch_size = batch_size_train,
  shuffle = TRUE
)

val_dl <- dataloader(
  val_dataset,
  batch_size = batch_size_val,
  shuffle = FALSE
)

message("Dataloaders created.")
message("Train batches: ", length(train_dl))
message("Val batches:   ", length(val_dl))

# metadata

input_dim   <- length(feature_cols)
num_species <- max(data_all_clean$label)  

message("Input dimension (num features): ", input_dim)
message("Number of species (output dim): ", num_species)

meta_list <- list(
  feature_cols = feature_cols,
  input_dim    = input_dim,
  num_species  = num_species
)

file_meta <- file.path(DIR_PROCESSED, "training_metadata_1000_species.rds")
saveRDS(meta_list, file_meta)
message("Saved training metadata to: ", file_meta)
