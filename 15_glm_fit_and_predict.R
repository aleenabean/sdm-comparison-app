# fit a single-species GLM SDM and project it to a lon/lat grid
# using the same env + coord encoding as the deep model

library(dplyr)
library(terra)

source("paths.R")
check_project_dirs()

deg2rad <- function(x) x * pi / 180

# fit GLM SDM
fit_glm_sdm <- function(df_glm,
                        response     = "presence",
                        exclude_cols = c("lon", "lat")) {
  
  if (!response %in% names(df_glm)) {
    stop("Response column '", response, "' not found in df_glm.")
  }
  
  predictors <- setdiff(names(df_glm), c(response, exclude_cols))
  if (length(predictors) == 0) {
    stop("No predictors left after excluding: ",
         paste(c(response, exclude_cols), collapse = ", "))
  }
  
  formula_str <- paste(response, "~", paste(predictors, collapse = " + "))
  form <- as.formula(formula_str)
  
  message("Fitting GLM with formula: ", formula_str)
  
  glm_fit <- glm(
    form,
    data   = df_glm,
    family = binomial()
  )
  
  glm_fit
}

# Predict GLM to a lon/lat grid & return a SpatRaster

# glm_model: object from fit_glm_sdm()
# env_stack: SpatRaster of WorldClim variables (same as used in 14)
# bbox: numeric vector or list: c(lon_min, lon_max, lat_min, lat_max)
# grid_res: grid resolution in degrees (e.g. 0.25)

# returns a list with:
# $raster: SpatRaster of probabilities
# $grid_df: data.frame with lon, lat, prob
#
predict_glm_to_grid <- function(glm_model,
                                env_stack,
                                bbox,
                                grid_res = 0.25) {

  if (is.numeric(bbox) && length(bbox) == 4) {
    lon_min <- bbox[1]
    lon_max <- bbox[2]
    lat_min <- bbox[3]
    lat_max <- bbox[4]
  } else if (is.list(bbox) &&
             all(c("lon_min", "lon_max", "lat_min", "lat_max") %in% names(bbox))) {
    lon_min <- bbox$lon_min
    lon_max <- bbox$lon_max
    lat_min <- bbox$lat_min
    lat_max <- bbox$lat_max
  } else {
    stop("bbox must be numeric c(lon_min, lon_max, lat_min, lat_max) or ",
         "a list with names lon_min, lon_max, lat_min, lat_max.")
  }
  
  lon_min <- max(lon_min, -180)
  lon_max <- min(lon_max,  180)
  lat_min <- max(lat_min,  -90)
  lat_max <- min(lat_max,   90)
  
  lon_seq <- seq(lon_min, lon_max, by = grid_res)
  lat_seq <- seq(lat_min, lat_max, by = grid_res)
  
  grid_df <- expand.grid(
    lon = lon_seq,
    lat = lat_seq
  )
  
  if (nrow(grid_df) == 0) {
    stop("Empty grid: check bbox and grid_res.")
  }
  
  message("Constructed grid with ", nrow(grid_df), " cells.")
  
  # extract env
  pts_grid <- terra::vect(
    grid_df,
    geom = c("lon", "lat"),
    crs  = "EPSG:4326"
  )
  
  env_vals <- terra::extract(env_stack, pts_grid, df = TRUE)
  env_vals <- env_vals[, -1, drop = FALSE]
  env_vals[is.na(env_vals)] <- 0
  
  # coord encoding
  lon_rad <- deg2rad(grid_df$lon)
  lat_rad <- deg2rad(grid_df$lat)
  
  newdata <- data.frame(
    lon = grid_df$lon,
    lat = grid_df$lat,
    x1  = sin(lon_rad),
    x2  = cos(lon_rad),
    x3  = sin(lat_rad),
    x4  = cos(lat_rad),
    env_vals
  )
  
  # match columns used by glm_model
  # extract predictor names from the model formula
  model_terms <- attr(terms(glm_model), "term.labels")
  
  # model_terms are like "x1", "x2", "wc2.1_2.5m_bio_1", ...
  missing_terms <- setdiff(model_terms, names(newdata))
  if (length(missing_terms) > 0) {
    stop("The following predictors used in the GLM are missing in newdata: ",
         paste(missing_terms, collapse = ", "))
  }
  
  # build the newdata for predict using exactly the columns the model expects
  newdata_for_predict <- newdata[, unique(c(model_terms)), drop = FALSE]
  
  prob_vec <- predict(glm_model, newdata = newdata_for_predict, type = "response")
  
  # build raster
  n_lon <- length(lon_seq)
  n_lat <- length(lat_seq)
  
  prob_mat <- matrix(NA_real_, nrow = n_lat, ncol = n_lon)
  
  lon_index <- match(grid_df$lon, lon_seq)
  lat_index <- match(grid_df$lat, lat_seq)
  
  for (i in seq_along(prob_vec)) {
    r <- lat_index[i]
    c <- lon_index[i]
    prob_mat[r, c] <- prob_vec[i]
  }
  
  rast_prob <- terra::rast(
    nrows = n_lat,
    ncols = n_lon,
    xmin  = lon_min,
    xmax  = lon_max,
    ymin  = lat_min,
    ymax  = lat_max,
    crs   = "EPSG:4326"
  )
  
  terra::values(rast_prob) <- as.vector(prob_mat)
  
  grid_out <- grid_df
  grid_out$prob <- prob_vec
  
  list(
    raster  = rast_prob,
    grid_df = grid_out
  )
}