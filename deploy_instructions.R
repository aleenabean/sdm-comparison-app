library(rsconnect)

# All SINR grid RDS files
sinr_grid_files <- Sys.glob("data_processed/sinr_grid_*.rds")

# Explicit list of files to bundle
app_files <- c(
  # Core app code
  "app.R",
  "paths.R",
  "14_glm_prepare_species_data.R",
  "15_glm_fit_and_predict.R",
  
  # Data actually used by the app
  "data_intermediate/inat_occurrences_1000_species.csv",
  "data_intermediate/iucn_ranges_app_species.rds",
  "data_intermediate/prototype_1000_species.csv",   # ← NEW: fix for your error
  
  "data_processed/training_points_1000_species.csv",
  "data_processed/sinr_1000species_iucn_eval.csv",
  
  # All precomputed SINR grids
  sinr_grid_files
)

# Deploy with *only* these files
rsconnect::deployApp(
  appDir   = ".",
  appFiles = app_files,
  appName  = "sdm_comparison_app",
  account  = "aleenamunshi",
  upload   = TRUE
)
