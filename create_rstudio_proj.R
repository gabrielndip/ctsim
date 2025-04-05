# create_rstudio_proj.R
# Setup script for ctsim package development environment

# Create a new RStudio project in the current directory
rstudioapi::initializeProject(path = ".")

# Initialize renv
renv::init()

# Install required packages for clinical trial simulation
renv::install(c(
  "tidyverse",    # Data manipulation 
  "survival",     # Survival analysis
  "flexsurv",     # Parametric survival models
  "brms",         # Bayesian regression models
  "msm",          # Multi-state models
  "testthat",     # Testing framework
  "styler",       # Code styling
  "lintr",        # Code linting
  "quarto",       # Documentation
  "here"          # Path management
))

# Snapshot dependencies
renv::snapshot()