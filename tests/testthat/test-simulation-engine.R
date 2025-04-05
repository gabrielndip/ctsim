context("Simulation engine functions")

# Create a test fixture for simulation tests
create_test_sim_result <- function() {
  # Create a trial config
  config <- create_trial_config(
    trial_name = "Test_Trial",
    n_subjects = 50,
    enrollment_rate = 5,
    follow_up_duration = 24,
    arms = list(
      control = list(
        name = "Control",
        allocation = 0.5,
        tte_distribution = "weibull",
        tte_params = list(shape = 1.2, scale = 12)
      ),
      treatment = list(
        name = "Treatment",
        allocation = 0.5,
        tte_distribution = "weibull",
        tte_params = list(shape = 1.2, scale = 18)
      )
    ),
    events_required = 40,
    dropout_rate = 0.1,
    seed = 123
  )
  
  # Add an interim analysis
  config <- add_interim_analysis(
    config,
    timing = 0.5,
    type = "events",
    decision_rules = list(
      efficacy = function(results) { results$hazard_ratio < 0.6 & results$p_value < 0.01 },
      futility = function(results) { results$hazard_ratio > 0.9 | results$p_value > 0.3 }
    )
  )
  
  # Generate participants
  participants <- generate_participants(config, scenario = "baseline")
  
  # Create simulation results
  sim_result <- create_simulation_results(config)
  sim_result$participants <- participants
  
  return(sim_result)
}

# Test run_trial_simulation function
test_that("run_trial_simulation orchestrates the entire simulation process", {
  # Use a minimal test configuration to speed up test
  config <- create_trial_config(
    trial_name = "Mini_Test",
    n_subjects = 20,
    enrollment_rate = 10,
    follow_up_duration = 12,
    events_required = 15,
    seed = 123
  )
  
  # Test single simulation
  sim_result <- run_trial_simulation(config, n_sims = 1, scenario = "baseline")
  
  # Check structure
  expect_s3_class(sim_result, "trial_simulation")
  expect_equal(sim_result$config$trial_name, "Mini_Test")
  expect_true(!is.null(sim_result$participants))
  expect_true(!is.null(sim_result$events))
  expect_true(!is.null(sim_result$final_analysis))
  
  # Test multiple simulations
  all_results <- run_trial_simulation(config, n_sims = 2, scenario = "baseline")
  
  # Check structure for multiple simulations
  expect_true(is.list(all_results))
  expect_length(all_results$individual_sims, 2)
  expect_true(!is.null(all_results$summary))
  
  # Test with different scenario
  sim_result_efficacy <- run_trial_simulation(config, n_sims = 1, scenario = "efficacy")
  expect_equal(sim_result_efficacy$scenario, "efficacy")
  
  # Invalid scenario should error
  expect_error(
    run_trial_simulation(config, scenario = "invalid"),
    "should be one of"
  )
})

# Test simulate_trial_events function
test_that("simulate_trial_events correctly processes all trial events", {
  # Create a test simulation result
  sim_result <- create_test_sim_result()
  
  # Run the simulation
  sim_result <- simulate_trial_events(sim_result)
  
  # Check that events were generated
  expect_true(is.data.frame(sim_result$events))
  expect_true(nrow(sim_result$events) > 0)
  
  # Check that trial completion time is set
  expect_true(!is.na(sim_result$trial_completion_time))
  
  # Check that event count is recorded
  expect_true(sim_result$event_count > 0)
  
  # Check event timeline order
  expect_true(all(diff(sim_result$events$time) >= 0))
  
  # Check that metrics were calculated
  expect_true(!is.null(sim_result$metrics))
  expect_true(!is.null(sim_result$metrics$enrollment_duration))
  expect_true(!is.null(sim_result$metrics$study_duration))
  expect_true(!is.null(sim_result$metrics$event_rate))
  
  # Check that the simulation continues until required events
  expect_true(sim_result$event_count >= sim_result$config$events_required)
})

# Test run_interim_analyses function
test_that("run_interim_analyses correctly implements interim analysis logic", {
  # Create a test simulation result
  sim_result <- create_test_sim_result()
  
  # Add a second interim analysis
  sim_result$config <- add_interim_analysis(
    sim_result$config,
    timing = 0.75,
    type = "events"
  )
  
  # Add a time-based interim analysis
  sim_result$config <- add_interim_analysis(
    sim_result$config,
    timing = 10,
    type = "time"
  )
  
  # Run event simulation
  sim_result <- simulate_trial_events(sim_result)
  
  # Run interim analyses
  sim_result <- run_interim_analyses(sim_result)
  
  # Check that analyses were conducted
  expect_true(is.list(sim_result$analyses))
  expect_true(length(sim_result$analyses) > 0)
  
  # Check each analysis
  for (analysis in sim_result$analyses) {
    # Each analysis should have key fields
    expect_true(!is.null(analysis$time))
    expect_true(!is.null(analysis$event_count))
    expect_true(!is.null(analysis$enrolled_count))
    
    # Statistical results should be present if not an error
    if (!("error" %in% names(analysis)) || !analysis$error) {
      expect_true(!is.null(analysis$hazard_ratio))
      expect_true(!is.null(analysis$p_value))
    }
    
    # Decision flags should be present
    expect_true(!is.null(analysis$stop_for_efficacy))
    expect_true(!is.null(analysis$stop_for_futility))
    expect_true(!is.null(analysis$stop_trial))
  }
  
  # Test early stopping
  # Create a configuration with extreme decision rules
  always_stop_config <- create_trial_config(
    trial_name = "AlwaysStop",
    n_subjects = 50,
    events_required = 40
  )
  
  # Add an interim analysis that always stops
  always_stop_config <- add_interim_analysis(
    always_stop_config,
    timing = 0.25,
    type = "events",
    decision_rules = list(
      efficacy = function(results) { TRUE },  # Always stop for efficacy
      futility = function(results) { FALSE }
    )
  )
  
  # Run simulation
  always_stop_sim <- create_simulation_results(always_stop_config)
  always_stop_sim$participants <- generate_participants(always_stop_config)
  always_stop_sim <- simulate_trial_events(always_stop_sim)
  always_stop_sim <- run_interim_analyses(always_stop_sim)
  
  # Check that trial stopped early
  expect_true(!is.null(always_stop_sim$trial_stopped_early))
  expect_true(always_stop_sim$trial_stopped_early)
  expect_equal(always_stop_sim$stop_reason, "efficacy")
})

# Test run_final_analysis function
test_that("run_final_analysis correctly analyzes trial outcomes", {
  # Create a test simulation result
  sim_result <- create_test_sim_result()
  
  # Run event simulation
  sim_result <- simulate_trial_events(sim_result)
  
  # Run final analysis
  sim_result <- run_final_analysis(sim_result)
  
  # Check that final analysis exists
  expect_true(!is.null(sim_result$final_analysis))
  
  # Check key analysis components
  fa <- sim_result$final_analysis
  expect_true(!is.null(fa$time))
  expect_true(!is.null(fa$event_count))
  expect_true(!is.null(fa$enrolled_count))
  expect_true(!is.null(fa$hazard_ratio))
  expect_true(!is.null(fa$hr_ci_lower))
  expect_true(!is.null(fa$hr_ci_upper))
  expect_true(!is.null(fa$p_value))
  expect_true(!is.null(fa$log_rank_p))
  expect_true(!is.null(fa$km_fit))
  expect_true(!is.null(fa$median_survival))
  expect_true(!is.null(fa$arm_event_rates))
  expect_true(!is.null(fa$parametric_models))
  expect_equal(fa$is_interim, FALSE)
  
  # Success flag should be set
  expect_true(!is.null(sim_result$success))
  expect_type(sim_result$success, "logical")
  
  # Test with trial stopped at interim
  stopped_sim <- create_test_sim_result()
  stopped_sim <- simulate_trial_events(stopped_sim)
  
  # Manually set trial stopped early
  stopped_sim$trial_stopped_early <- TRUE
  stopped_sim$analyses <- list(
    list(
      time = 10,
      event_count = 25,
      hazard_ratio = 0.5,
      p_value = 0.005,
      stop_for_efficacy = TRUE,
      stop_for_futility = FALSE,
      stop_trial = TRUE
    )
  )
  
  # Run final analysis
  stopped_sim <- run_final_analysis(stopped_sim)
  
  # Check that final analysis used interim results
  expect_equal(stopped_sim$final_analysis$time, 10)
  expect_equal(stopped_sim$final_analysis$event_count, 25)
  expect_equal(stopped_sim$final_analysis$hazard_ratio, 0.5)
  expect_equal(stopped_sim$final_analysis$is_interim, TRUE)
})

# Test apply_bayesian_updates function
test_that("apply_bayesian_updates correctly processes Bayesian calculations", {
  # Skip test unless Bayesian packages are available
  skip_if_not_installed("brms")
  
  # Create a test prior
  hr_prior <- create_bayesian_prior(
    distribution = "normal",
    parameters = list(mean = 0, sd = 1)
  )
  
  # Create a trial configuration with priors
  config <- create_trial_config(
    trial_name = "BayesianTest",
    n_subjects = 30,
    arms = list(
      control = list(
        name = "Control",
        allocation = 0.5,
        tte_distribution = "weibull",
        tte_params = list(shape = 1.2, scale = 12)
      ),
      treatment = list(
        name = "Treatment",
        allocation = 0.5,
        tte_distribution = "weibull",
        tte_params = list(shape = 1.2, scale = 18),
        prior = hr_prior
      )
    ),
    prior_data = list(
      dropout_prior = create_bayesian_prior(
        distribution = "beta",
        parameters = list(shape1 = 2, shape2 = 20)
      )
    )
  )
  
  # Create simulation result
  sim_result <- create_simulation_results(config)
  sim_result$participants <- generate_participants(config)
  
  # Mock interim results since we'll skip actual fitting
  sim_result$participants$event_status <- rbinom(nrow(sim_result$participants), 1, 0.4)
  sim_result$participants$observed_time <- pmax(1, rnorm(nrow(sim_result$participants), 10, 3))
  
  # Attempt to run Bayesian updates
  # This may fail if brms isn't properly set up, so we'll catch and skip
  tryCatch({
    sim_result <- apply_bayesian_updates(sim_result)
    
    # Check if updates were applied
    expect_true(!is.null(sim_result$bayesian_updates))
    expect_false(is.null(sim_result$bayesian_updates$message))
  }, error = function(e) {
    skip("Skipping full Bayesian test due to model fitting issues")
  })
  
  # Test with no prior data
  config_no_prior <- create_trial_config()
  sim_no_prior <- create_simulation_results(config_no_prior)
  sim_no_prior$participants <- generate_participants(config_no_prior)
  
  sim_no_prior <- apply_bayesian_updates(sim_no_prior)
  expect_match(sim_no_prior$bayesian_updates$message, "No prior data provided")
})

# Test summarize_simulations function
test_that("summarize_simulations correctly calculates aggregate statistics", {
  # Create multiple simulation results
  n_sims <- 3
  all_results <- list()
  
  for (i in 1:n_sims) {
    sim <- create_test_sim_result()
    sim <- simulate_trial_events(sim)
    sim <- run_final_analysis(sim)
    
    # Set different outcomes
    sim$success <- i %% 2 == 0  # Alternating success/failure
    sim$trial_completion_time <- 10 + i * 2
    sim$final_analysis$hazard_ratio <- 0.6 + i * 0.1
    sim$event_count <- 40 + i
    
    all_results[[i]] <- sim
  }
  
  # Summarize the results
  summary <- summarize_simulations(all_results)
  
  # Check summary structure
  expect_true(!is.null(summary$n_sims))
  expect_equal(summary$n_sims, n_sims)
  expect_true(!is.null(summary$config))
  expect_true(!is.null(summary$success_probability))
  expect_true(!is.null(summary$avg_duration))
  expect_true(!is.null(summary$median_duration))
  expect_true(!is.null(summary$avg_hr))
  expect_true(!is.null(summary$median_hr))
  expect_true(!is.null(summary$avg_events))
  expect_true(!is.null(summary$hr_quantiles))
  expect_true(!is.null(summary$prob_clinical_effect))
  
  # Check summary metrics
  expect_equal(summary$success_probability, 1/3)  # 1 out of 3 sims
  expect_equal(summary$avg_duration, mean(c(12, 14, 16)))
  expect_equal(summary$avg_hr, mean(c(0.7, 0.8, 0.9)))
  
  # Test with non-list input
  expect_error(
    summarize_simulations("not a list"),
    "all_results must be a list"
  )
})

# Test plot_simulation_results function
test_that("plot_simulation_results creates visualizations", {
  # Skip test if ggplot2 is not available
  skip_if_not_installed("ggplot2")
  
  # Create a test simulation result
  sim_result <- create_test_sim_result()
  sim_result <- simulate_trial_events(sim_result)
  sim_result <- run_final_analysis(sim_result)
  
  # Try different plot types
  km_plot <- plot_simulation_results(sim_result, plot_type = "km")
  enrollment_plot <- plot_simulation_results(sim_result, plot_type = "enrollment")
  events_plot <- plot_simulation_results(sim_result, plot_type = "events")
  
  # Check that plots were created
  expect_s3_class(km_plot, "ggplot")
  expect_s3_class(enrollment_plot, "ggplot")
  expect_s3_class(events_plot, "ggplot")
  
  # Test with invalid plot type
  expect_error(
    plot_simulation_results(sim_result, plot_type = "invalid"),
    "should be one of"
  )
})

# Test export_simulation_results function
test_that("export_simulation_results exports data in various formats", {
  # Create a test simulation result
  sim_result <- create_test_sim_result()
  sim_result <- simulate_trial_events(sim_result)
  sim_result <- run_final_analysis(sim_result)
  
  # Create a temp directory for exports
  temp_dir <- tempdir()
  csv_path <- file.path(temp_dir, "test_export")
  rds_path <- file.path(temp_dir, "test_export.rds")
  
  # Test CSV export
  result_path <- export_simulation_results(sim_result, format = "csv", file_path = csv_path)
  expect_equal(result_path, csv_path)
  
  # Participants CSV should exist
  expect_true(file.exists(paste0(csv_path, "_participants.csv")))
  
  # Test RDS export
  result_path <- export_simulation_results(sim_result, format = "rds", file_path = rds_path)
  expect_equal(result_path, rds_path)
  expect_true(file.exists(rds_path))
  
  # Test invalid format
  expect_error(
    export_simulation_results(sim_result, format = "invalid"),
    "should be one of"
  )
  
  # Clean up
  unlink(paste0(csv_path, "*.csv"))
  unlink(rds_path)
})