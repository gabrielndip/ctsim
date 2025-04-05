context("Visualization functions")

# Create test fixture for visualization tests
create_test_sim_for_vis <- function() {
  # Create a trial config
  config <- create_trial_config(
    trial_name = "VisTest",
    n_subjects = 40,
    enrollment_rate = 5,
    follow_up_duration = 12,
    arms = list(
      control = list(
        name = "Control",
        allocation = 0.5,
        tte_distribution = "weibull",
        tte_params = list(shape = 1.2, scale = 6)
      ),
      treatment = list(
        name = "Treatment",
        allocation = 0.5,
        tte_distribution = "weibull",
        tte_params = list(shape = 1.2, scale = 9)
      )
    ),
    events_required = 25,
    dropout_rate = 0.05,
    seed = 789
  )
  
  # Generate participants
  participants <- generate_participants(config, scenario = "baseline")
  
  # Create simulation results
  sim_result <- create_simulation_results(config)
  sim_result$participants <- participants
  
  # Run simulation steps
  sim_result <- simulate_trial_events(sim_result)
  sim_result <- run_final_analysis(sim_result)
  
  return(sim_result)
}

# Test plot_km_curve function
test_that("plot_km_curve creates Kaplan-Meier plots", {
  # Skip if ggplot2 or gridExtra are not available
  skip_if_not_installed(c("ggplot2", "gridExtra"))
  
  # Create test data
  sim_result <- create_test_sim_for_vis()
  
  # Create plots with different options
  plot1 <- plot_km_curve(sim_result, conf_int = TRUE, risk_table = FALSE)
  plot2 <- plot_km_curve(sim_result, conf_int = FALSE, risk_table = FALSE)
  
  # With risk table, it returns a grid.arrange object
  plot3 <- plot_km_curve(sim_result, conf_int = TRUE, risk_table = TRUE)
  
  # Check plot types
  expect_s3_class(plot1, "ggplot")
  expect_s3_class(plot2, "ggplot")
  
  # plot3 is a grid.arrange when risk_table=TRUE
  # but this is hard to test due to grid's internal structure
  # so we'll skip detailed checks on that
})

# Test plot_enrollment function
test_that("plot_enrollment creates enrollment plots", {
  # Skip if ggplot2 is not available
  skip_if_not_installed("ggplot2")
  
  # Create test data
  sim_result <- create_test_sim_for_vis()
  
  # Create plots with different options
  plot1 <- plot_enrollment(sim_result, by_arm = TRUE)
  plot2 <- plot_enrollment(sim_result, by_arm = FALSE)
  
  # Check plot types
  expect_s3_class(plot1, "ggplot")
  expect_s3_class(plot2, "ggplot")
  
  # Check titles
  expect_equal(plot1$labels$title, "Enrollment Over Time")
  
  # Custom title
  plot3 <- plot_enrollment(sim_result, title = "My Custom Title")
  expect_equal(plot3$labels$title, "My Custom Title")
})

# Test plot_events function
test_that("plot_events creates event plots", {
  # Skip if ggplot2 is not available
  skip_if_not_installed("ggplot2")
  
  # Create test data
  sim_result <- create_test_sim_for_vis()
  
  # Create plots with different options
  plot1 <- plot_events(sim_result, by_arm = TRUE)
  plot2 <- plot_events(sim_result, by_arm = FALSE)
  plot3 <- plot_events(sim_result, include_censoring = TRUE)
  
  # Check plot types
  expect_s3_class(plot1, "ggplot")
  expect_s3_class(plot2, "ggplot")
  expect_s3_class(plot3, "ggplot")
  
  # Check titles
  expect_equal(plot1$labels$title, "Events Over Time")
  
  # Custom title
  plot4 <- plot_events(sim_result, title = "Event Analysis")
  expect_equal(plot4$labels$title, "Event Analysis")
})

# Test plot_bayesian_posterior function
test_that("plot_bayesian_posterior creates posterior plots", {
  # Skip if ggplot2 is not available
  skip_if_not_installed("ggplot2")
  
  # Create mock Bayesian posterior
  mock_posterior <- list(
    samples = list(
      b_armtreatment = rnorm(1000, -0.2, 0.1)  # Log hazard ratio
    ),
    median = -0.2,
    mean = -0.19
  )
  
  # Create simulation result with mock Bayesian updates
  sim_result <- create_test_sim_for_vis()
  sim_result$bayesian_updates <- list(
    hr_posterior = mock_posterior,
    prob_efficacy = 0.85,
    prob_superiority = 0.95
  )
  
  # Create hazard ratio plot
  hr_plot <- plot_bayesian_posterior(sim_result, parameter = "hr")
  
  # Create efficacy probability plot
  efficacy_plot <- plot_bayesian_posterior(sim_result, parameter = "efficacy")
  
  # Check plot types
  expect_s3_class(hr_plot, "ggplot")
  expect_s3_class(efficacy_plot, "ggplot")
  
  # Default title for HR plot
  expect_equal(hr_plot$labels$title, "Posterior Distribution of Hazard Ratio")
  
  # Custom title
  custom_plot <- plot_bayesian_posterior(sim_result, parameter = "hr", 
                                       title = "Treatment Effect Distribution")
  expect_equal(custom_plot$labels$title, "Treatment Effect Distribution")
  
  # Test with missing Bayesian updates
  sim_result_no_bayes <- create_test_sim_for_vis()
  sim_result_no_bayes$bayesian_updates <- NULL
  
  fallback_plot <- plot_bayesian_posterior(sim_result_no_bayes, parameter = "hr")
  # Should return a placeholder plot, still a ggplot object
  expect_s3_class(fallback_plot, "ggplot")
})

# Test plot_simulation_summary function
test_that("plot_simulation_summary creates summary plots for multiple simulations", {
  # Skip if ggplot2 is not available
  skip_if_not_installed("ggplot2")
  
  # Create mock simulation summary
  sim_summary <- list(
    n_sims = 100,
    success_probability = 0.82,
    early_stopping_probability = 0.35,
    study_durations = rnorm(100, 24, 3),
    hazard_ratios = rlnorm(100, log(0.75), 0.2)
  )
  
  # Create different summary plots
  success_plot <- plot_simulation_summary(sim_summary, metric = "success_rate")
  stopping_plot <- plot_simulation_summary(sim_summary, metric = "early_stopping")
  duration_plot <- plot_simulation_summary(sim_summary, metric = "study_duration")
  hr_plot <- plot_simulation_summary(sim_summary, metric = "hazard_ratio")
  
  # Check plot types
  expect_s3_class(success_plot, "ggplot")
  expect_s3_class(stopping_plot, "ggplot")
  expect_s3_class(duration_plot, "ggplot")
  expect_s3_class(hr_plot, "ggplot")
  
  # Check default titles
  expect_equal(success_plot$labels$title, "Probability of Trial Success")
  expect_equal(hr_plot$labels$title, "Distribution of Estimated Hazard Ratios")
  
  # Custom title
  custom_plot <- plot_simulation_summary(sim_summary, metric = "success_rate", 
                                      title = "Success Probability")
  expect_equal(custom_plot$labels$title, "Success Probability")
  
  # Test with missing data
  sim_summary_missing <- list(success_probability = 0.75)
  limited_plot <- plot_simulation_summary(sim_summary_missing, metric = "study_duration")
  # Should return a placeholder plot, still a ggplot object
  expect_s3_class(limited_plot, "ggplot")
})

# Test create_simulation_dashboard function
test_that("create_simulation_dashboard combines multiple plots", {
  # Skip if required packages are not available
  skip_if_not_installed(c("ggplot2", "gridExtra", "grid"))
  
  # Create test data
  sim_result <- create_test_sim_for_vis()
  
  # Create dashboard with default plots
  dashboard <- create_simulation_dashboard(sim_result)
  
  # Create dashboard with specific plots
  selected_dashboard <- create_simulation_dashboard(
    sim_result,
    plots = c("km", "enrollment"),
    ncol = 1
  )
  
  # It's difficult to test grid.arrange objects directly
  # but we can check they're created without error
  expect_true(!is.null(dashboard))
  expect_true(!is.null(selected_dashboard))
})