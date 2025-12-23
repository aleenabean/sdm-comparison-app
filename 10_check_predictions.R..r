#inspect predictions of the SINR-style model for a few sample points

library(dplyr)
library(readr)
library(torch)

source("paths.R")
check_project_dirs()

# load metadata, data, species lookup

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

file_processed <- file.path(DIR_PROCESSED, "training_points_1000_species.csv")
if (!file.exists(file_processed)) {
  stop("Processed training file not found: ", file_processed,
       "\nRun 06_build_training_table.R first.")
}
data_all <- read_csv(file_processed, show_col_types = FALSE)
message("Loaded processed data: ", nrow(data_all), " rows.")

file_prototype <- file.path(DIR_INTERMEDIATE, "prototype_1000_species.csv")
if (!file.exists(file_prototype)) {
  stop("Prototype species list not found: ", file_prototype,
       "\nRun 04_select_1000_species.R first.")
}
prototype_species <- read_csv(file_prototype, show_col_types = FALSE)

species_lookup <- prototype_species %>%
  mutate(label = prototype_index + 1L) %>%
  select(label, iucn_scientific_name = scientific_name, taxon_group)

# load model

source("08_model_and_loss.R")  # defines sinr_mlp

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

# sample some train points with complete features

set.seed(123)
n_examples <- 10

sample_rows <- data_all %>%
  filter(split == "train") %>%
  filter(complete.cases(across(all_of(feature_cols)))) %>%
  slice_sample(n = n_examples)

missing_feats <- setdiff(feature_cols, names(sample_rows))
if (length(missing_feats) > 0) {
  stop("Missing feature columns in sample_rows: ",
       paste(missing_feats, collapse = ", "))
}
if (!("prototype_index" %in% names(sample_rows))) {
  stop("prototype_index column missing in sample_rows.")
}

sample_rows <- sample_rows %>%
  mutate(label = prototype_index + 1L)

# run model

feat_mat <- as.matrix(sample_rows[, feature_cols, drop = FALSE])
x <- torch_tensor(feat_mat, dtype = torch_float(), device = device)

p <- model(x)  # [n_examples, num_species], already probabilities (sigmoid)

p_cpu <- p$to(device = "cpu")
p_mat <- as.matrix(as_array(p_cpu))

#print top-5 predictions for each example

top_k <- 5

for (i in seq_len(n_examples)) {
  row_i <- sample_rows[i, ]
  true_label <- as.integer(row_i$label)
  
  true_sp <- species_lookup %>%
    filter(label == true_label) %>%
    slice(1)
  
  true_name  <- true_sp$iucn_scientific_name
  true_group <- true_sp$taxon_group
  
  cat("\n----------------------------------------------------\n")
  cat("Example", i, "\n")
  cat("  Location: lon =", round(row_i$decimalLongitude, 4),
      ", lat =", round(row_i$decimalLatitude, 4), "\n")
  cat("  True species:", true_name, "(", true_group, ")\n")
  
  probs_i <- p_mat[i, ]
  
  ord <- order(probs_i, decreasing = TRUE)
  top_idx <- ord[1:top_k]
  top_probs <- probs_i[top_idx]
  
  pred_tbl <- tibble(
    label = top_idx,
    prob  = top_probs
  ) %>%
    left_join(species_lookup, by = "label") %>%
    arrange(desc(prob))
  
  cat("  Top", top_k, "predictions:\n")
  for (j in seq_len(nrow(pred_tbl))) {
    cat(sprintf("    #%d: %s (%s)  p = %.4f\n",
                j,
                pred_tbl$iucn_scientific_name[j],
                pred_tbl$taxon_group[j],
                pred_tbl$prob[j]))
  }
}
cat("\nDone.\n")
