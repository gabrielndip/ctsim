# Function to create project directory structure
create_project_structure <- function(base_dir = ".") {
  # Main directories
  dirs <- c(
    "R", "data/inputs", "data/sample",
    "tests/testthat", "vignettes", "reports"
  )

  # Create directories
  for (dir in dirs) {
    dir_path <- file.path(base_dir, dir)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
      message("Created directory: ", dir_path)
    }
  }

  # Create R script files
  r_files <- c(
    "survival_models.R", "bayesian_utils.R",
    "simulation_engine.R", "data_generation.R", "visualization.R"
  )

  for (file in r_files) {
    file_path <- file.path(base_dir, "R", file)
    if (!file.exists(file_path)) {
      file.create(file_path)
      # Add standard header
      writeLines(
        c(
          "# Clinical Trial Simulation",
          paste0("# ", file),
          "# Created: ", format(Sys.Date(), "%Y-%m-%d"),
          "# Description: [Add description]",
          "",
          "library(tidyverse)",
          ""
        ),
        file_path
      )
      message("Created file: ", file_path)
    }
  }

  # Create config files
  writeLines(
    "linters: with_defaults(\n  line_length_linter(120),\n  object_length_linter(40),\n  object_name_linter(styles = c(\"snake_case\")),\n  assignment_linter = NULL\n)",
    file.path(base_dir, ".lintr")
  )

  writeLines(
    "style_pkg(scope = \"spaces\", strict = TRUE, indent_by = 2)",
    file.path(base_dir, ".styler")
  )

  writeLines(
    c(".Rproj.user", ".Rhistory", ".RData", ".Ruserdata", "renv/library", "renv/local", "renv/cellar", "renv/staging", ".DS_Store"),
    file.path(base_dir, ".gitignore")
  )

  # Create README
  readme_content <- c(
    "# Clinical Trial Simulation",
    "",
    "## Overview",
    "A statistical simulation framework for clinical trials focusing on time-to-event analysis, Bayesian updates, and participant retention modeling.",
    "",
    "## Features",
    "- Parametric survival models for time-to-event simulation",
    "- Bayesian updating for incorporating new clinical data",
    "- Dropout prediction models for trial retention",
    "- Multi-state modeling for participant journey",
    "- Scenario comparison for trial planning",
    "",
    "## Getting Started",
    "1. Clone this repository",
    "2. Open the R project",
    "3. Run `renv::restore()` to install dependencies",
    "4. See vignettes/ directory for examples"
  )

  writeLines(readme_content, file.path(base_dir, "README.md"))
  message("Created README.md")

  # Create DESCRIPTION file
  desc_content <- c(
    "Package: ctsim",
    "Title: Clinical Trial Simulation Framework",
    "Version: 0.0.1",
    "Authors@R: person('Gabriel Teku', email = 'gabbyteku@gmail.com', role = c('aut', 'cre'))",
    "Description: A framework for simulating clinical trials with survival analysis and Bayesian updating.",
    "License: MIT",
    "Encoding: UTF-8",
    "LazyData: true",
    "Roxygen: list(markdown = TRUE)",
    "RoxygenNote: 7.2.3",
    "Imports:",
    "    tidyverse,",
    "    survival,",
    "    flexsurv,",
    "    brms,",
    "    msm,",
    "    here"
  )

  writeLines(desc_content, file.path(base_dir, "DESCRIPTION"))
  message("Created DESCRIPTION file")
}

# Run the function
create_project_structure()
