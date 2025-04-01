# Create a new RStudio project in the current directory
rstudioapi::initializeProject(path = ".")

# Initialize renv
renv::init()

# Install packages
renv::install(c(
  "tidyverse", "survival", "flexsurv", "brms", 
  "msm", "testthat", "styler", "lintr", "quarto", "here"
))

# Snapshot dependencies
renv::snapshot()
