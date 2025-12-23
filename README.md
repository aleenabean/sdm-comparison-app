Comparing Single-Species and Multi-Species Species Distribution Models

This repository contains the code for an interactive Shiny application that compares a traditional single-species Species Distribution Model (SDM) with a multi-species neural SDM inspired by Spatial Implicit Neural Representations (SINR).

Live app:
https://aleenamunshi.shinyapps.io/sdm_comparison_app/

PROJECT OVERVIEW:

Species distribution models predict where species are likely to occur based on environmental covariates and observed occurrences. Classical SDMs typically fit one model per species, while newer machine-learning approaches train multi-species models that share information across taxa.
This project compares:
* a single-species GLM SDM fit independently for each species using pseudo-absences, and
* a multi-species neural model trained offline and evaluated via precomputed probability grids.

Both models are evaluated on the same spatial grid and compared against IUCN Red List expert range polygons using AUC and average precision (AP).

The project emphasizes:
* methodological comparison between single- and multi-species SDMs,
* shared evaluation against expert range maps
* and practical deployment under real-world computational constraints.


REPOSITORY STRUCTURE:

This repository intentionally includes only the files needed to run the app and inspect the modeling logic.

Files included:
- app.R
- Main Shiny application.
- paths.R (Centralized path and directory management)
- 14_glm_prepare_species_data.R (Data preparation for the single-species GLM.)
- 15_glm_fit_and_predict.R (GLM fitting and grid-based prediction.)
- deploy_instructions.R (Notes related to deployment and app structure.)
- .gitignore
- .rscignore

WHAT IS NOT INCLUDED:

The following were run offline and are intentionally excluded:
* raw data ingestion and cleaning scripts,
* large intermediate datasets,
* neural network training code,
* full global raster datasets,
* precomputation scripts for SINR grids.

These steps are described conceptually in the accompanying write-up and reflected in the app’s behavior, but excluded to keep the repository lightweight and reproducible.

MODELS:

Single-Species GLM:

- Logistic regression (binomial GLM)
- Presence-only iNaturalist/GBIF data
- Pseudo-absences sampled from buffered IUCN bounding boxes*
- WorldClim 2.1 bioclimatic variables
- Evaluated on a fixed 0.25-degree grid

Multi-Species SINR-Style Model:

- Multi-layer perceptron with shared parameters across species
- Continuous spatial encoding (latitude, longitude, environment)
- Single-positive multi-label loss
- Trained offline
- Probability grids precomputed per species
- Loaded and cropped on demand in the Shiny app

EVALUATION:

Both models are evaluated against IUCN expert range polygons:
- Grid cells inside the polygon are labeled positive
- Cells outside (within the evaluation region) are labeled negative
  
Metrics:
- AUC (Area Under the ROC Curve)
- Average Precision (AP)
Evaluation is range-based rather than true presence–absence, and applied consistently across models.

RUNNING THE APP LOCALLY:

Requirements:
- R (>= 4.2 recommended)
- Packages: shiny, terra, sf, dplyr, ggplot2, pROC, viridis, bslib

Steps:

1. Clone the repository
git clone https://github.com/aleenabean/sdm-comparison-app.git
cd sdm-comparison-app

2. Ensure required data files are available (paths managed in paths.R)

3. Run the app
shiny::runApp("app.R")

WorldClim data are downloaded lazily when models are run to avoid large startup overhead.

REFERENCES:

Cole, E., et al. (2023).
Spatial Implicit Neural Representations for Global Species Mapping.
arXiv:2306.02564

Harrell, L., et al. (2025).
A Heterogeneous Graph Neural Network for Species Distribution Modeling.
arXiv:2503.11900
