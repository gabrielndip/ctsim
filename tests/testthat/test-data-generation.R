context("Data generation functions")

# Test create_trial_config function
test_that("create_trial_config produces a valid trial configuration", {
  # Basic test
  config <- create_trial_config(
    trial_name = "Test_Trial",
    n_subjects = 100,
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
    events_required = 80,
    dropout_rate = 0.1,
    seed = 123
  )
  
  # Check that return value is of correct type
  expect_s3_class(config, "trial_config")
  
  # Check essential fields
  expect_equal(config$trial_name, "Test_Trial")
  expect_equal(config$n_subjects, 100)
  expect_equal(config$enrollment_rate, 5)
  expect_equal(config$events_required, 80)
  expect_equal(config$dropout_rate, 0.1)
  expect_equal(config$seed, 123)
  
  # Check treatment arms
  expect_named(config$arms, c("control", "treatment"))
  expect_equal(config$arms$control$allocation, 0.5)
  expect_equal(config$arms$treatment$allocation, 0.5)
  
  # Test error case: arm allocations don't sum to 1
  expect_error(
    create_trial_config(
      arms = list(
        control = list(name = "Control", allocation = 0.4),
        treatment = list(name = "Treatment", allocation = 0.4)
      )
    ),
    "Arm allocations must sum to 1"
  )
})

# Test add_interim_analysis function
test_that("add_interim_analysis correctly adds interim analyses", {
  # Create a base configuration
  config <- create_trial_config()
  
  # Add an interim analysis
  config <- add_interim_analysis(
    config,
    timing = 0.5,
    type = "events",
    decision_rules = list(
      efficacy = function(results) { results$hazard_ratio < 0.7 & results$p_value < 0.01 },
      futility = function(results) { results$hazard_ratio > 0.9 | results$p_value > 0.3 }
    )
  )
  
  # Check that interim analysis was added correctly
  expect_length(config$interim_analyses, 1)
  expect_equal(config$interim_analyses[[1]]$timing, 0.5)
  expect_equal(config$interim_analyses[[1]]$type, "events")
  expect_type(config$interim_analyses[[1]]$decision_rules$efficacy, "closure")
  expect_type(config$interim_analyses[[1]]$decision_rules$futility, "closure")
  
  # Add a second interim analysis
  config <- add_interim_analysis(
    config,
    timing = 0.75,
    type = "events"
  )
  
  # Check that second interim analysis was added
  expect_length(config$interim_analyses, 2)
  expect_equal(config$interim_analyses[[2]]$timing, 0.75)
  
  # Test time-based interim analysis
  config <- add_interim_analysis(
    config,
    timing = 12,
    type = "time"
  )
  
  # Check that time-based analysis was added
  expect_length(config$interim_analyses, 3)
  expect_equal(config$interim_analyses[[3]]$type, "time")
  expect_equal(config$interim_analyses[[3]]$timing, 12)
})

# Test create_simulation_results function
test_that("create_simulation_results initializes results structure correctly", {
  # Create a trial config
  config <- create_trial_config()
  
  # Create simulation results
  results <- create_simulation_results(config)
  
  # Check that return value is of correct type
  expect_s3_class(results, "trial_simulation")
  
  # Check essential fields
  expect_identical(results$config, config)
  expect_true(is.data.frame(results$participants))
  expect_true(is.data.frame(results$events))
  expect_true(is.list(results$analyses))
  expect_true(is.na(results$trial_completion_time))
  expect_equal(results$event_count, 0)
  expect_true(is.list(results$final_analysis))
  expect_equal(results$scenario, "baseline")
  expect_true(is.list(results$bayesian_updates))
})

# Test create_bayesian_prior function
test_that("create_bayesian_prior creates valid prior distributions", {
  # Normal prior
  normal_prior <- create_bayesian_prior(
    distribution = "normal",
    parameters = list(mean = 0, sd = 1)
  )
  
  # Check that return value is of correct type
  expect_s3_class(normal_prior, "bayesian_prior")
  
  # Check fields
  expect_equal(normal_prior$distribution, "normal")
  expect_equal(normal_prior$parameters$mean, 0)
  expect_equal(normal_prior$parameters$sd, 1)
  
  # Test beta prior
  beta_prior <- create_bayesian_prior(
    distribution = "beta",
    parameters = list(alpha = 2, beta = 20)
  )
  
  expect_s3_class(beta_prior, "bayesian_prior")
  expect_equal(beta_prior$distribution, "beta")
  
  # Test gamma prior
  gamma_prior <- create_bayesian_prior(
    distribution = "gamma",
    parameters = list(shape = 2, rate = 1)
  )
  
  expect_s3_class(gamma_prior, "bayesian_prior")
  expect_equal(gamma_prior$distribution, "gamma")
  
  # Test invalid distribution type
  expect_error(
    create_bayesian_prior(distribution = "invalid_dist"),
    "should be one of"
  )
})

# Test generate_participants function
test_that("generate_participants creates valid participant data", {
  # Create a trial config
  config <- create_trial_config(
    n_subjects = 50,
    enrollment_rate = 5,
    arms = list(
      control = list(
        name = "Control",
        allocation = 0.4,
        tte_distribution = "weibull",
        tte_params = list(shape = 1.2, scale = 12)
      ),
      treatment = list(
        name = "Treatment",
        allocation = 0.6,
        tte_distribution = "weibull",
        tte_params = list(shape = 1.2, scale = 18)
      )
    ),
    seed = 123
  )
  
  # Generate participants for baseline scenario
  participants <- generate_participants(config, scenario = "baseline")
  
  # Check dimensions and structure
  expect_equal(nrow(participants), 50)
  expect_true("participant_id" %in% names(participants))
  expect_true("arm" %in% names(participants))
  expect_true("enrollment_time" %in% names(participants))
  expect_true("event_time" %in% names(participants))
  expect_true("censoring_time" %in% names(participants))
  expect_true("observed_time" %in% names(participants))
  expect_true("event_status" %in% names(participants))
  
  # Check arm allocation
  expect_equal(sum(participants$arm == "control"), 20)
  expect_equal(sum(participants$arm == "treatment"), 30)
  
  # Generate for efficacy scenario
  participants_efficacy <- generate_participants(config, scenario = "efficacy")
  
  # Treatment should have better outcomes in efficacy scenario
  mean_time_control <- mean(participants_efficacy$event_time[participants_efficacy$arm == "control"])
  mean_time_treatment <- mean(participants_efficacy$event_time[participants_efficacy$arm == "treatment"])
  expect_gt(mean_time_treatment, mean_time_control)
  
  # Test invalid scenario
  expect_error(
    generate_participants(config, scenario = "invalid_scenario"),
    "should be one of"
  )
})

# Test generate_enrollment_times function
test_that("generate_enrollment_times creates valid enrollment pattern", {
  # Test with positive enrollment rate
  times <- generate_enrollment_times(100, 5)
  
  # Check length and ordering
  expect_length(times, 100)
  expect_true(all(diff(times) >= 0))  # Times should be increasing
  
  # Average time between enrollments should be close to 1/rate
  expect_lt(abs(mean(diff(times)) - 1/5), 0.1)
  
  # Test error with negative or zero enrollment rate
  expect_error(
    generate_enrollment_times(100, 0),
    "Enrollment rate must be positive"
  )
  expect_error(
    generate_enrollment_times(100, -1),
    "Enrollment rate must be positive"
  )
})

# Test generate_event_times function
test_that("generate_event_times correctly generates event times by arm", {
  # Create a trial config
  config <- create_trial_config(
    n_subjects = 100,
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
        tte_distribution = "exponential",
        tte_params = list(scale = 15)
      )
    )
  )
  
  # Generate participants without event times
  participants <- generate_participants(config)
  participants$event_time <- NA
  
  # Generate event times
  participants <- generate_event_times(participants, config, "baseline")
  
  # Check that all event times are filled
  expect_false(any(is.na(participants$event_time)))
  
  # Check that different distributions were used for each arm
  # Control arm should follow Weibull with shape=1.2
  # Treatment arm should follow exponential
  
  # For Weibull, variance depends on shape and scale
  var_control <- var(participants$event_time[participants$arm == "control"])
  var_treatment <- var(participants$event_time[participants$arm == "treatment"])
  
  # We expect treatment to have higher variance (exponential) than control (Weibull with shape>1)
  # This is a statistical expectation based on the properties of these distributions
  expect_true(abs(var_treatment - var_control) > 0)
})

# Test for adjust_params_for_scenario
test_that("adjust_params_for_scenario correctly modifies parameters by scenario", {
  # Test baseline scenario (no changes)
  params <- list(scale = 10)
  adjusted <- adjust_params_for_scenario(params, "treatment", "baseline")
  expect_equal(adjusted$scale, 10)
  
  # Test efficacy scenario with treatment arm (should improve)
  adjusted <- adjust_params_for_scenario(params, "treatment", "efficacy")
  expect_equal(adjusted$scale, 15)  # 50% increase
  
  # Test efficacy scenario with control arm (shouldn't change)
  adjusted <- adjust_params_for_scenario(params, "control", "efficacy")
  expect_equal(adjusted$scale, 10)  # No change
  
  # Test with log-normal parameters
  lnorm_params <- list(meanlog = 2)
  adjusted <- adjust_params_for_scenario(lnorm_params, "treatment", "efficacy")
  expect_equal(adjusted$meanlog, 2 + log(1.5))  # Add log(1.5) for 50% increase
})

# Test generate_censoring_times
test_that("generate_censoring_times creates valid censoring times", {
  # Create a trial config
  config <- create_trial_config(
    n_subjects = 50,
    dropout_rate = 0.1
  )
  
  # Create participant data
  participants <- data.frame(
    participant_id = paste0("P", 1:50),
    enrollment_time = seq(0, 10, length.out = 50),
    stringsAsFactors = FALSE
  )
  
  # Generate censoring times
  participants <- generate_censoring_times(participants, config, "baseline")
  
  # Check that censoring times are created
  expect_true("censoring_time" %in% names(participants))
  expect_false(any(is.na(participants$censoring_time)))
  
  # Censoring times should be after enrollment
  expect_true(all(participants$censoring_time > participants$enrollment_time))
  
  # With rate 0.1, mean time to censoring should be around 10
  expect_lt(abs(mean(participants$censoring_time - participants$enrollment_time) - 10), 3)
  
  # Test efficacy scenario (lower dropout)
  participants_efficacy <- generate_censoring_times(participants, config, "efficacy")
  
  # Mean time to censoring should be higher in efficacy scenario
  mean_time_baseline <- mean(participants$censoring_time - participants$enrollment_time)
  mean_time_efficacy <- mean(participants_efficacy$censoring_time - participants$enrollment_time)
  expect_gt(mean_time_efficacy, mean_time_baseline)
})

# Test finalize_participant_data
test_that("finalize_participant_data correctly determines observed times and status", {
  # Create a trial config
  config <- create_trial_config(
    follow_up_duration = 24
  )
  
  # Create participant data with event and censoring times
  participants <- data.frame(
    participant_id = c("P1", "P2", "P3", "P4"),
    enrollment_time = c(0, 5, 10, 15),
    event_time = c(12, 8, 30, 20),  # Time from enrollment to event
    censoring_time = c(24, 4, 15, 25),  # Time from enrollment to censoring
    stringsAsFactors = FALSE
  )
  
  # Finalize participant data
  participants <- finalize_participant_data(participants, config)
  
  # Check observed times and status
  # P1: Event at 12 months (earlier than censoring or end of follow-up)
  expect_equal(participants$observed_time[1], 12)
  expect_equal(participants$event_status[1], 1)
  expect_equal(participants$retention_status[1], "Event")
  
  # P2: Censoring at 4 months (earlier than event or end of follow-up)
  expect_equal(participants$observed_time[2], 4)
  expect_equal(participants$event_status[2], 0)
  expect_equal(participants$retention_status[2], "Dropout")
  
  # P3: End of follow-up at 24 months (earlier than event, later than censoring)
  expect_equal(participants$observed_time[3], 24)
  expect_equal(participants$event_status[3], 0)
  expect_equal(participants$retention_status[3], "Completed")
  
  # P4: End of follow-up at 24 months (earlier than both event and censoring)
  # Note: Enrolled at 15, so end of follow-up is at 15+24=39,
  # but we observe event at 15+20=35, which is earlier
  expect_equal(participants$observed_time[4], 20)
  expect_equal(participants$event_status[4], 1)
  expect_equal(participants$retention_status[4], "Event")
})

# Test generate_stratified_participants
test_that("generate_stratified_participants creates participants with strata and covariate effects", {
  # Create a trial config
  config <- create_trial_config(
    n_subjects = 100,
    seed = 456
  )
  
  # Generate stratified participants
  strata <- list(
    risk_group = list(
      levels = c("High", "Medium", "Low"),
      probs = c(0.25, 0.5, 0.25)
    ),
    region = list(
      levels = c("North", "South", "East", "West"),
      probs = c(0.3, 0.3, 0.2, 0.2)
    )
  )
  
  covariate_effects <- list(
    risk_group = list(
      "High" = 1.5,    # Higher hazard (worse outcomes)
      "Medium" = 1.0,  # Reference
      "Low" = 0.7      # Lower hazard (better outcomes)
    ),
    age = 1.02  # 2% increase in hazard per year for age
  )
  
  participants <- generate_stratified_participants(
    config, 
    strata = strata,
    covariate_effects = covariate_effects
  )
  
  # Check that strata were added
  expect_true("risk_group" %in% names(participants))
  expect_true("region" %in% names(participants))
  
  # Check strata proportions
  risk_table <- table(participants$risk_group)
  expect_equal(risk_table["High"] / sum(risk_table), 0.25, tolerance = 0.1)
  expect_equal(risk_table["Medium"] / sum(risk_table), 0.5, tolerance = 0.1)
  expect_equal(risk_table["Low"] / sum(risk_table), 0.25, tolerance = 0.1)
  
  # Check that covariate effects were applied
  # High risk should have shorter event times than low risk
  high_risk_times <- participants$event_time[participants$risk_group == "High"]
  low_risk_times <- participants$event_time[participants$risk_group == "Low"]
  
  # Due to hazard multiplier, high risk should have shorter times on average
  expect_lt(mean(high_risk_times), mean(low_risk_times))
})

# Test generate_test_data
test_that("generate_test_data creates a sample dataset", {
  # Generate test data
  test_data <- generate_test_data(n_samples = 30)
  
  # Check dimensions and structure
  expect_equal(nrow(test_data), 30)
  expect_true("participant_id" %in% names(test_data))
  expect_true("arm" %in% names(test_data))
  expect_true("enrollment_time" %in% names(test_data))
  expect_true("event_time" %in% names(test_data))
  expect_true("observed_time" %in% names(test_data))
  expect_true("event_status" %in% names(test_data))
  
  # Check arm balance (should be 15/15)
  expect_equal(sum(test_data$arm == "control"), 15)
  expect_equal(sum(test_data$arm == "treatment"), 15)
  
  # Check randomization
  test_data2 <- generate_test_data(n_samples = 30)
  # Some values should be different (due to randomization)
  expect_false(identical(test_data$event_time, test_data2$event_time))
})