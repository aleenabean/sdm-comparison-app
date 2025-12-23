# FINAL!!!!
# Shiny app to compare single-species GLM SDM vs multi-species SINR model

library(shiny)
library(bslib)
library(dplyr)
library(readr)
library(terra)
library(sf)
library(pROC)
library(ggplot2)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)
library(googledrive)

source("paths.R")
source("14_glm_prepare_species_data.R")
source("15_glm_fit_and_predict.R")
check_project_dirs()

deg2rad <- function(x) x * pi / 180

# fixed grid resolution for BOTH models
GRID_RES_GLOBAL <- 0.25

if (!exists("WORLDCLIM_DRIVE_ID")) {
  WORLDCLIM_DRIVE_ID <- "1iYwuMdbGTSMgPuY0VYMPBZq6kP2TPSLk"
}

# download + cache WorldClim in tempdir() for this R session and build a SpatRaster
load_env_stack_for_app <- function() {
  wc_dir <- file.path(tempdir(), "worldclim_bio_app")
  if (!dir.exists(wc_dir)) {
    dir.create(wc_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  tif_files <- list.files(
    wc_dir,
    pattern   = "\\.tif$",
    full.names = TRUE,
    recursive  = TRUE
  )
  
  if (length(tif_files) == 0) {
    message("Downloading WorldClim data from Google Drive (will take a couple minutes)...")
    cat("Downloading WorldClim data from Google Drive (will take a couple minutes)...\n")
    
    tmp_zip <- file.path(tempdir(), "worldclim_wc2_1_2_5m_bio.zip")
    
    googledrive::drive_deauth()
    googledrive::drive_download(
      googledrive::as_id(WORLDCLIM_DRIVE_ID),
      path      = tmp_zip,
      overwrite = TRUE
    )
    
    utils::unzip(tmp_zip, exdir = wc_dir)
    unlink(tmp_zip)
    
    tif_files <- list.files(
      wc_dir,
      pattern   = "\\.tif$",
      full.names = TRUE,
      recursive  = TRUE
    )
    
    if (length(tif_files) == 0) {
      stop("WorldClim download finished but no .tif files were found in ", wc_dir)
    }
  } else {
    message("Using cached WorldClim data from tempdir().")
  }
  
  terra::rast(tif_files)
}


# IUCN ranges (small app-specific RDS)
FILE_IUCN_RDS <- file.path(DIR_INTERMEDIATE, "iucn_ranges_app_species.rds")

if (!file.exists(FILE_IUCN_RDS)) {
  stop(
    "IUCN app-species RDS not found at: ", FILE_IUCN_RDS,
    "\nRun 20_make_small_iucn.R in development to create it."
  )
}

iucn_all <- readRDS(FILE_IUCN_RDS)
if (!inherits(iucn_all, "sf")) stop("IUCN object is not sf.")
if (!"sci_name" %in% names(iucn_all)) stop("IUCN sf lacks 'sci_name' column.")

# prototype species (for SINR labels)
FILE_PROTO <- file.path(DIR_INTERMEDIATE, "prototype_1000_species.csv")
prototype_species <- readr::read_csv(FILE_PROTO, show_col_types = FALSE)

#species list

inat_occ <- load_inat_occurrences_default()

sp_col_inat <- "iucn_scientific_name"
if (!sp_col_inat %in% names(inat_occ)) {
  stop(
    "Column 'iucn_scientific_name' not found in inat_occ.\n",
    "Got columns: ", paste(names(inat_occ), collapse = ", ")
  )
}

species_choices <- {
  inat_species   <- unique(na.omit(inat_occ[[sp_col_inat]]))
  proto_species  <- prototype_species$scientific_name
  iucn_species   <- unique(iucn_all$sci_name)
  intersect(intersect(inat_species, proto_species), iucn_species) |>
    sort()
}
if (length(species_choices) == 0) {
  stop("No overlapping species between iNat, prototypes, and IUCN.")
}

#cross-species SINR eval

FILE_SINR_EVAL <- file.path(DIR_PROCESSED, "sinr_1000species_iucn_eval.csv")
sinr_eval <- if (file.exists(FILE_SINR_EVAL)) {
  readr::read_csv(FILE_SINR_EVAL, show_col_types = FALSE)
} else {
  NULL
}

if (!is.null(sinr_eval)) {
  obs_counts <- inat_occ %>%
    group_by(.data[[sp_col_inat]]) %>%
    summarise(n_obs = n(), .groups = "drop")
  names(obs_counts)[1] <- "scientific_name"
  
  sinr_eval <- sinr_eval %>%
    left_join(obs_counts, by = "scientific_name")
}

max_pos_global <- if (!is.null(sinr_eval) && "n_pos_all" %in% names(sinr_eval)) {
  max(sinr_eval$n_pos_all, na.rm = TRUE)
} else 1000

max_obs_global <- if (!is.null(sinr_eval) && "n_obs" %in% names(sinr_eval)) {
  max(sinr_eval$n_obs, na.rm = TRUE)
} else 1000
if (!is.finite(max_obs_global) || max_obs_global <= 0) max_obs_global <- 1000

# ------------------------------------------------------------
# WORLD BASEMAP FOR RANGE SUMMARY
# ------------------------------------------------------------

world_basemap <- rnaturalearth::ne_countries(
  scale       = "medium",
  returnclass = "sf"
)
if (sf::st_crs(world_basemap)$epsg != 4326) {
  world_basemap <- sf::st_transform(world_basemap, 4326)
}

#helping functions

get_label_index <- function(sp_name) {
  row <- prototype_species %>%
    filter(scientific_name == sp_name) %>%
    slice(1)
  if (nrow(row) == 0) return(NA_integer_)
  row$prototype_index + 1L
}

get_iucn_polygon <- function(sp_name) {
  if (is.null(iucn_all)) {
    return(NULL)
  }
  
  iucn_sp_sf <- iucn_all[iucn_all$sci_name == sp_name, ]
  if (nrow(iucn_sp_sf) == 0) return(NULL)
  
  iucn_sp_sf <- sf::st_make_valid(iucn_sp_sf)
  iucn_sp_sf <- iucn_sp_sf %>% dplyr::summarise()
  
  if (is.na(sf::st_crs(iucn_sp_sf))) {
    sf::st_crs(iucn_sp_sf) <- 4326
  } else if (!is.null(sf::st_crs(iucn_sp_sf)$epsg) &&
             sf::st_crs(iucn_sp_sf)$epsg != 4326) {
    iucn_sp_sf <- sf::st_transform(iucn_sp_sf, 4326)
  }
  
  iucn_sp_sf
}

get_species_region_info <- function(sp_name) {
  iucn_sp_sf <- get_iucn_polygon(sp_name)
  if (is.null(iucn_sp_sf)) {
    return(list(
      continents = NA_character_,
      countries  = NA_character_
    ))
  }
  
  if (sf::st_crs(iucn_sp_sf)$epsg != 4326) {
    iucn_sp_sf <- sf::st_transform(iucn_sp_sf, 4326)
  }
  
  idx_list <- sf::st_intersects(iucn_sp_sf, world_basemap, sparse = TRUE)
  idx <- idx_list[[1]]
  
  if (length(idx) == 0) {
    return(list(
      continents = NA_character_,
      countries  = NA_character_
    ))
  }
  
  countries  <- sort(unique(world_basemap$admin[idx]))
  continents <- sort(unique(world_basemap$continent[idx]))
  
  list(
    continents = continents,
    countries  = countries
  )
}

eval_grid_against_iucn <- function(sp_name, grid_df, prob_col = "prob") {
  iucn_sp_sf <- get_iucn_polygon(sp_name)
  if (is.null(iucn_sp_sf)) {
    return(list(
      auc   = NA_real_,
      ap    = NA_real_,
      n_pos = NA_integer_,
      n_neg = NA_integer_
    ))
  }
  
  if (!all(c("lon", "lat", prob_col) %in% names(grid_df))) {
    stop("grid_df must have columns lon, lat, and ", prob_col)
  }
  
  pts_sf <- sf::st_as_sf(grid_df, coords = c("lon", "lat"), crs = 4326)
  
  inside_list <- sf::st_intersects(pts_sf, iucn_sp_sf, sparse = TRUE)
  y_true <- ifelse(lengths(inside_list) > 0, 1L, 0L)
  prob   <- grid_df[[prob_col]]
  
  n_pos <- sum(y_true == 1L)
  n_neg <- sum(y_true == 0L)
  
  if (length(unique(y_true)) < 2 || n_pos == 0L || n_neg == 0L) {
    return(list(
      auc   = NA_real_,
      ap    = NA_real_,
      n_pos = n_pos,
      n_neg = n_neg
    ))
  }
  
  roc_obj <- pROC::roc(response = y_true, predictor = prob, quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  ord <- order(prob, decreasing = TRUE)
  y_sorted <- y_true[ord]
  cum_tp   <- cumsum(y_sorted == 1L)
  positions <- seq_along(y_sorted)
  precisions <- cum_tp / positions
  
  num_pos <- sum(y_sorted == 1L)
  if (num_pos == 0L) {
    ap_val <- NA_real_
  } else {
    ap_val <- sum(precisions[y_sorted == 1L]) / num_pos
  }
  
  list(
    auc   = auc_val,
    ap    = ap_val,
    n_pos = n_pos,
    n_neg = n_neg
  )
}

# SINR predictions from precomputed grids (so no torch yay!)

slugify_species <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

predict_sinr_to_grid <- function(sp_name, bbox) {
  slug <- slugify_species(sp_name)
  file_rds <- file.path(DIR_PROCESSED, paste0("sinr_grid_", slug, ".rds"))
  
  if (!file.exists(file_rds)) {
    stop(
      "Precomputed SINR grid not found for species: ", sp_name, "\n",
      "Expected file: ", file_rds
    )
  }
  
  obj <- readRDS(file_rds)
  if (is.null(obj$grid_df)) {
    stop(
      "Precomputed SINR file for ", sp_name,
      " does not contain a 'grid_df' component."
    )
  }
  
  grid_df <- obj$grid_df
  
  lon_min <- bbox[1]; lon_max <- bbox[2]
  lat_min <- bbox[3]; lat_max <- bbox[4]
  
  inside <- grid_df$lon >= lon_min & grid_df$lon <= lon_max &
    grid_df$lat >= lat_min & grid_df$lat <= lat_max
  
  grid_df_sub <- grid_df[inside, , drop = FALSE]
  
  if (nrow(grid_df_sub) == 0) {
    grid_df_sub <- grid_df
  }
  
  list(
    grid_df = grid_df_sub
  )
}

# UI

theme_lora <- bs_theme(
  version    = 5,
  bootswatch = "darkly",
  base_font  = font_google("Lora")
)

ui <- fluidPage(
  theme = theme_lora,
  titlePanel("Single-species GLM vs Multi-species SINR SDM"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "species",
        "Choose species:",
        choices  = species_choices,
        selected = species_choices[1]
      ),
      numericInput(
        "n_pseudo_abs",
        "GLM pseudo-absences (max):",
        value = 1000,
        min   = 100,
        max   = 5000,
        step  = 100
      ),
      actionButton("run_button", "Run both models"),
      br(),
      helpText(
        paste0(
          "Map grid resolution is fixed at ",
          GRID_RES_GLOBAL,
          "° (same for GLM and SINR)."
        )
      ),
      helpText("AUC/AP are computed vs IUCN polygons on this grid.")
    ),
    mainPanel(
      tabsetPanel(
        id = "main_tabs",
        
        tabPanel(
          "Maps & Models",
          h4("Maps"),
          helpText("Note: the black outline on each map shows the IUCN expert range polygon."),
          fluidRow(
            column(6, plotOutput("glmPlot",  height = "350px")),
            column(6, plotOutput("sinrPlot", height = "350px"))
          ),
          hr(),
          h4("Metrics vs IUCN (grid-based)"),
          tableOutput("metricsTable"),
          hr(),
          h4("Range Summary"),
          verbatimTextOutput("speciesSummary"),
          hr(),
          h4("Programmer tools"),
          helpText("Cross-species SINR summary, simple inference benchmark, and CSV download."),
          conditionalPanel(
            condition = "true",
            sliderInput(
              "min_pos_eval",
              "Minimum positive grid cells (n_pos_all) to include:",
              min   = 0,
              max   = max_pos_global,
              value = min(20, max_pos_global),
              step  = 10
            ),
            sliderInput(
              "max_obs",
              "Maximum iNat observations (n_obs) to include (focus on rarer species):",
              min   = 0,
              max   = max_obs_global,
              value = max_obs_global,
              step  = 50
            ),
            tableOutput("crossSpeciesTable"),
            plotOutput("crossSpeciesPlot", height = "300px")
          ),
          hr(),
          h4("Model efficiency (inference only, current species)"),
          actionButton("bench_button", "Benchmark GLM vs SINR"),
          verbatimTextOutput("benchOutput"),
          hr(),
          h4("Download SINR per-species metrics"),
          downloadButton("download_sinr_eval", "Download CSV")
        ),
        
        tabPanel(
          "Background & Context",
          h3("What this app is doing"),
          p("This app compares a traditional single-species Species Distribution Model (SDM) ",
            "against a multi-species deep-learning SDM inspired by Spatial Implicit Neural ",
            "Representations (SINR)."),
          p("Both models take presence-only occurrence data plus environmental covariates ",
            "and try to learn where a species is likely to be present across space. ",
            "Predictions are evaluated on a fixed 0.25° grid, using IUCN expert polygons ",
            "as a proxy for true range."),
          h3("Why SDMs are hard"),
          tags$ul(
            tags$li("Most biodiversity data are presence-only: we know where a species ",
                    "was observed, not where it is truly absent."),
            tags$li("Sampling is strongly biased toward certain regions, taxa, and ",
                    "accessible locations. Models can accidentally learn sampling bias, ",
                    "not biological reality."),
            tags$li("Evaluation is tricky: high-quality presence–absence data are rare, ",
                    "so researchers often rely on expert maps or pseudo-absences.")
          ),
          h3("Single-species vs multi-species SDMs"),
          p("Traditional SDMs fit one model per species, ignoring information about other ",
            "species. Multi-species approaches instead learn a shared representation of ",
            "environmental space and can 'borrow strength' across species."),
          p("This can especially help rare species with very few observations, because ",
            "their model can share parameters with better-sampled species that live in ",
            "similar environments."),
          h3("How we measure performance (AUC and AP)"),
          tags$ul(
            tags$li(strong("AUC (area under ROC curve): "),
                    "probability that the model ranks a random presence higher ",
                    "than a random absence. 0.5 ≈ random, 1.0 = perfect."),
            tags$li(strong("Average Precision (AP): "),
                    "summarizes precision across probability thresholds, with more weight ",
                    "on the highest-probability predictions, which is especially relevant ",
                    "for rare species.")
          ),
          p("In this app, AUC and AP are computed on a grid of prediction cells, using ",
            "IUCN polygons as the 'positive' class and cells outside those polygons as ",
            "negatives within the evaluation region.")
        ),
        
        tabPanel(
          "Research (GNN & SINR)",
          h3("SINR-inspired multi-species model (this app)"),
          p("This app implements a small prototype inspired by the Spatial Implicit Neural ",
            "Representations (SINR) model of Cole et al. It is designed for demonstration ",
            "and teaching, not for production conservation decisions."),
          h4("How our SINR-style model works"),
          tags$ul(
            tags$li("Input: a 4D trigonometric encoding of latitude/longitude plus 19 WorldClim ",
                    "Bioclim variables (and elevation if available)."),
            tags$li("Architecture: multi-layer perceptron with residual connections, ",
                    "outputting probabilities for 1000 species simultaneously."),
            tags$li("Loss function: single-positive multi-label loss. For each occurrence, ",
                    "the observed species is treated as the only known positive; all other ",
                    "species at that location and random background locations are treated as ",
                    "pseudo-negatives with a higher weight on the positive term.")
          ),
          h4("How this differs from the full SINR model"),
          tags$ul(
            tags$li("Scale: this prototype uses ~1000 species and tens of thousands of observations; ",
                    "Cole et al. train on ~47,000 species and 35M+ iNaturalist observations globally."),
            tags$li("Computation: this app trains on a single laptop (Apple M3) with relatively few ",
                    "epochs and minimal hyperparameter tuning; the original work uses large-scale GPU training."),
            tags$li("Loss and regularization: the loss here is a simplified variant with a single ",
                    "positive-weight parameter; the original paper carefully tunes multiple terms, regularization, ",
                    "and sampling strategies."),
            tags$li("Coverage: this prototype is restricted to mammals and amphibians with IUCN polygons; ",
                    "the full SINR model covers many more taxa and regions.")
          ),
          h4("Prototype disclaimer"),
          p("The predictions and metrics shown here are for methodological comparison only. ",
            "They inherit biases from presence-only data, pseudo-absence sampling, and IUCN polygons, ",
            "and should not be used directly for conservation or management decisions."),
          hr(),
          h3("Heterogeneous GNN comparison"),
          p("Harrell et al. (2025) propose a heterogeneous graph neural network (GNN) for species ",
            "distribution modeling. Their model builds a graph with species nodes and location nodes, ",
            "and learns to predict detection edges between them."),
          h4("Their data (NCEAS benchmark) vs ours"),
          tags$ul(
            tags$li("GNN: trained and evaluated on the NCEAS benchmark dataset, which provides ",
                    "presence-only training data and independent presence–absence survey data ",
                    "for regions such as Ontario, South America, Switzerland, New South Wales, ",
                    "Australian Wet Tropics, and New Zealand."),
            tags$li("This app: uses GBIF/iNaturalist presence-only observations, WorldClim covariates, ",
                    "and IUCN Red List polygons as a range-based evaluation target.")
          ),
          h4("Model differences"),
          tags$ul(
            tags$li("GNN: explicitly models a species–location graph and uses message passing between ",
                    "species and locations to share information."),
            tags$li("This app: uses an MLP with an implicit continuous representation of space, without ",
                    "explicit graph structure."),
            tags$li("GNN evaluation: AUCROC against held-out presence–absence surveys in each region."),
            tags$li("This app: AUC and AP against IUCN polygons on a coarse grid.")
          ),
          p("Because the data sources, objectives, and evaluation targets differ, results in this app ",
            "are not directly comparable to the GNN results. Instead, the app is meant to illustrate ",
            "the general idea of multi-species models (like SINR and the GNN) versus traditional ",
            "single-species SDMs.")
        ),
        
        tabPanel(
          "Credits & Contact",
          h3("Data sources"),
          tags$ul(
            tags$li("IUCN Red List: expert range polygons for mammals and amphibians (2025-2)."),
            tags$li("GBIF / iNaturalist: presence-only occurrence records used for training."),
            tags$li("WorldClim 2.1: 19 Bioclim variables (and elevation) at 2.5-arcminute resolution.")
          ),
          h3("Methods & references"),
          tags$ul(
            tags$li("Cole et al. (2023) – Spatial Implicit Neural Representations (SINR) for global species mapping."),
            tags$li("Harrell et al. (2025) – heterogeneous GNN for species distribution modeling on the NCEAS benchmark."),
            tags$li("Standard single-species SDM practices (GLM with pseudo-absences) as a baseline.")
          ),
          h3("Authorship"),
          p("Shiny app and SINR-style prototype implementation by ", strong("Aleena Munshi"), "."),
          p("This is a ", strong("Stats 20 Final Project"), "."),
          p("Source code available on GitHub: ",tags$a(href = "https://github.com/aleenabean/sdm-comparison-app",
              target = "_blank","github.com/aleenabean/sdm-comparison-app")),
          h3("Contact"),
          p("Email: ", tags$a(href = "mailto:aleenamunshi001@ucla.edu", "aleenamunshi001@ucla.edu"))
        )
      )
    )
  )
)

# server

server <- function(input, output, session) {
  
  # cache for WorldClim SpatRaster (download once per R session)
  env_stack_rv <- reactiveVal(NULL)
  
  # run both models for the chosen species
  run_results <- eventReactive(input$run_button, {
    sp_name <- input$species
    n_bg    <- as.integer(input$n_pseudo_abs)
    
    #  load WorldClim when the user actually runs models
    if (is.null(env_stack_rv())) {
      withProgress(message = "Downloading WorldClim data (will take a couple minutes)...", value = 0, {
        env_stack_rv(load_env_stack_for_app())
      })
    }
    env_stack <- env_stack_rv()
    
    # bbox from IUCN polygon
    iucn_sp_sf <- get_iucn_polygon(sp_name)
    if (!is.null(iucn_sp_sf)) {
      bb <- sf::st_bbox(iucn_sp_sf)
      margin <- 1
      bbox <- c(
        max(as.numeric(bb["xmin"]) - margin, -180),
        min(as.numeric(bb["xmax"]) + margin,  180),
        max(as.numeric(bb["ymin"]) - margin,  -90),
        min(as.numeric(bb["ymax"]) + margin,   90)
      )
    } else {
      bbox <- c(-20, 20, -20, 20)
    }
    
    region_info <- get_species_region_info(sp_name)
    
    # iNat obs used for training for this species
    n_obs_species <- inat_occ %>%
      filter(.data[[sp_col_inat]] == sp_name) %>%
      nrow()
    
    # GLM data prep, fit, predict, eval
    df_glm <- prepare_glm_data_for_species(
      scientific_name = sp_name,
      n_pseudo_abs    = n_bg,
      env_stack       = env_stack,
      inat_occ        = inat_occ,
      bbox_margin     = 1
    )
    glm_fit <- fit_glm_sdm(df_glm)
    
    proj_glm <- predict_glm_to_grid(
      glm_model = glm_fit,
      env_stack = env_stack,
      bbox      = bbox,
      grid_res  = GRID_RES_GLOBAL
    )
    eval_glm <- eval_grid_against_iucn(
      sp_name,
      proj_glm$grid_df,
      prob_col = "prob"
    )
    glm_res <- list(
      raster = proj_glm$raster,
      eval   = eval_glm,
      model  = glm_fit
    )
    
    # SINR predict and eval on same bbox
    proj_sinr <- predict_sinr_to_grid(sp_name, bbox)
    
    eval_sinr <- eval_grid_against_iucn(
      sp_name,
      proj_sinr$grid_df,
      prob_col = "prob"
    )
    
    sinr_res <- list(
      grid_df = proj_sinr$grid_df,
      eval    = eval_sinr
    )
    
    list(
      species = sp_name,
      glm     = glm_res,
      sinr    = sinr_res,
      bbox    = bbox,
      region  = region_info,
      n_obs   = n_obs_species
    )
  })
  
  #maps
  
  output$glmPlot <- renderPlot({
    res <- run_results()
    if (is.null(res) || is.null(res$glm)) return(NULL)
    
    r       <- res$glm$raster
    sp_name <- res$species
    iucn_sf <- get_iucn_polygon(sp_name)
    
    if (!is.null(iucn_sf)) {
      iucn_sv <- terra::vect(iucn_sf)
      e <- terra::ext(iucn_sv)
      margin <- 1
      e_buf <- e + c(-margin, margin, -margin, margin)
      
      r_crop <- terra::crop(r, e_buf)
      r_plot <- terra::disagg(r_crop, fact = 2, method = "bilinear")
      
      terra::plot(
        r_plot,
        col  = viridis::viridis(200),
        main = "GLM Probability of Presence",
        xlab = "Longitude (°)",
        ylab = "Latitude (°)",
        axes = TRUE
      )
      plot(iucn_sv, add = TRUE, border = "black", lwd = 0.8)
    } else {
      r_plot <- terra::disagg(r, fact = 2, method = "bilinear")
      terra::plot(
        r_plot,
        col  = viridis::viridis(200),
        main = "GLM Probability of Presence",
        xlab = "Longitude (°)",
        ylab = "Latitude (°)",
        axes = TRUE
      )
    }
  })
  
  output$sinrPlot <- renderPlot({
    res <- run_results()
    if (is.null(res) || is.null(res$sinr)) return(NULL)
    
    grid_df <- res$sinr$grid_df
    if (is.null(grid_df) || nrow(grid_df) == 0) {
      plot.new()
      title("SINR map: no predictions available")
      return(invisible())
    }
    
    sp_name <- res$species
    
    iucn_sf <- NULL
    try({
      iucn_sf <- get_iucn_polygon(sp_name)
    }, silent = TRUE)
    
    p <- ggplot(grid_df, aes(x = lon, y = lat, fill = prob)) +
      geom_raster() +
      scale_fill_viridis_c(option = "viridis") +
      labs(
        title = "SINR Probability of Presence",
        x     = "Longitude (°)",
        y     = "Latitude (°)",
        fill  = "Prob."
      ) +
      theme_minimal()
    
    if (!is.null(iucn_sf)) {
      p <- p +
        geom_sf(
          data        = iucn_sf,
          inherit.aes = FALSE,
          fill        = NA,
          color       = "black",
          linewidth   = 0.4
        ) +
        coord_sf()
    } else {
      p <- p + coord_equal()
    }
    
    print(p)
  })
  
  #metrics & range summary
  
  output$metricsTable <- renderTable({
    res <- run_results()
    if (is.null(res)) return(NULL)
    
    data.frame(
      Model = c("GLM", "SINR (multi-species)"),
      AUC   = round(c(res$glm$eval$auc,  res$sinr$eval$auc), 3),
      AP    = round(c(res$glm$eval$ap,   res$sinr$eval$ap),  3),
      stringsAsFactors = FALSE
    )
  })
  
  output$speciesSummary <- renderText({
    res <- run_results()
    if (is.null(res)) return("Run models to see range summary.")
    
    cont <- res$region$continents
    ctr  <- res$region$countries
    
    cont_str <- if (is.null(cont) || all(is.na(cont))) {
      "Not summarized by continent in this prototype (range shown on maps)."
    } else {
      paste(cont, collapse = ", ")
    }
    
    ctr_str <- if (is.null(ctr) || all(is.na(ctr))) {
      "Not summarized by country in this prototype."
    } else if (length(ctr) > 10) {
      paste0(paste(ctr[1:10], collapse = ", "),
             ", ... (", length(ctr), " countries total)")
    } else {
      paste(ctr, collapse = ", ")
    }
    
    paste0(
      "Continents: ", cont_str, "\n",
      "Countries:  ", ctr_str
    )
  })
  
  # programmer tools (cross-species summary)
  
  filtered_sinr <- reactive({
    if (is.null(sinr_eval)) return(NULL)
    df <- sinr_eval
    if (!is.null(input$min_pos_eval)) {
      df <- df %>% filter(is.na(n_pos_all) | n_pos_all >= input$min_pos_eval)
    }
    if (!is.null(input$max_obs) && "n_obs" %in% names(df)) {
      df <- df %>% filter(is.na(n_obs) | n_obs <= input$max_obs)
    }
    df
  })
  
  output$crossSpeciesTable <- renderTable({
    df <- filtered_sinr()
    if (is.null(df)) {
      return(data.frame(
        message = "sinr_1000species_iucn_eval.csv not found in data_processed/.",
        stringsAsFactors = FALSE
      ))
    }
    df %>%
      select(
        scientific_name, n_obs, n_pos_all, n_neg_all, auc, ap, note
      ) %>%
      arrange(desc(auc)) %>%
      head(50)
  })
  
  output$crossSpeciesPlot <- renderPlot({
    df <- filtered_sinr()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    ggplot(df, aes(x = log10(pmax(n_obs, 1)), y = auc)) +
      geom_point(alpha = 0.5) +
      labs(
        x = "log10(iNat observations + 1)",
        y = "SINR AUC vs IUCN",
        title = "SINR performance across species"
      ) +
      theme_minimal()
  })
  
  #inference benchmark
  
  bench_results <- eventReactive(input$bench_button, {
    res <- run_results()
    if (is.null(res)) {
      return("Please run the models for a species first.")
    }
    
    #making sure env is loaded
    if (is.null(env_stack_rv())) {
      env_stack_rv(load_env_stack_for_app())
    }
    env_stack <- env_stack_rv()
    
    sp_name <- res$species
    bbox    <- res$bbox
    
    # time GLM prediction
    t_glm <- system.time({
      predict_glm_to_grid(
        glm_model = res$glm$model,
        env_stack = env_stack,
        bbox      = bbox,
        grid_res  = GRID_RES_GLOBAL
      )
    })["elapsed"]
    
    # time SINR prediction
    t_sinr <- system.time({
      predict_sinr_to_grid(
        sp_name,
        bbox = bbox
      )
    })["elapsed"]
    
    paste0(
      "Inference time on same grid for ", sp_name, ":\n",
      "  GLM:  ", round(as.numeric(t_glm), 3),  " seconds\n",
      "  SINR: ", round(as.numeric(t_sinr), 3), " seconds"
    )
  })
  
  output$benchOutput <- renderText({
    bench_results()
  })
  
  #download SINR metrics
  
  output$download_sinr_eval <- downloadHandler(
    filename = function() {
      "sinr_1000species_iucn_eval.csv"
    },
    content = function(file) {
      if (is.null(sinr_eval)) {
        stop("sinr_1000species_iucn_eval.csv not available in data_processed/.")
      }
      readr::write_csv(sinr_eval, file)
    }
  )
}

shinyApp(ui = ui, server = server)

