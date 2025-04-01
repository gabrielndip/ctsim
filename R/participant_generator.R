# Synthetic Participant Data Generator
# Functions to generate synthetic participant data for clinical trial simulations

#' @title Generate synthetic participant data for a clinical trial
#' 
#' @param config Trial configuration object
#' @param scenario Scenario type ("baseline", "efficacy", "update")
#' @return Data frame of participant data
generate_participants <- function(
    config,
    scenario = c("baseline", "efficacy", "update")
) {
  scenario <- match.arg(scenario)
  
  # Set random seed for reproducibility
  set.seed(config$seed)
  
  # Calculate number of participants per arm
  arm_counts <- sapply(config$arms, function(arm) {
    round(arm$allocation * config$n_subjects)
  })
  
  # Adjust if rounding causes total to be different from n_subjects
  diff <- config$n_subjects - sum(arm_counts)
  if (diff != 0) {
    # Add/subtract the difference to/from the largest arm
    largest_arm <- which.max(arm_counts)
    arm_counts[largest_arm] <- arm_counts[largest_arm] + diff
  }
  
  # Generate participant IDs
  participant_ids <- paste0("P", sprintf("%04d", 1:config$n_subjects))
  
  # Generate arm assignments
  arm_names <- names(config$arms)
  arm_assignments <- rep(arm_names, arm_counts)
  
  # Basic demographic information (can be expanded)
  age <- round(rnorm(config$n_subjects, mean = 55, sd = 12))
  age <- ifelse(age < 18, 18, ifelse(age > 85, 85, age))  # Clip to realistic range
  
  gender <- sample(c("Male", "Female"), config$n_subjects, replace = TRUE, prob = c(0.48, 0.52))
  
  # Create base participant data frame
  participants <- data.frame(
    participant_id = participant_ids,
    arm = arm_assignments,
    age = age,
    gender = gender,
    enrollment_time = NA,
    event_time = NA,
    censoring_time = NA,
    observed_time = NA,
    event_status = NA,
    retention_status = NA,
    stringsAsFactors = FALSE
  )
  
  # Generate enrollment times (staggered enrollment)
  enrollment_rate <- config$enrollment_rate
  participants$enrollment_time <- generate_enrollment_times(config$n_subjects, enrollment_rate)
  
  # Generate event times based on arm-specific distributions
  participants <- generate_event_times(participants, config, scenario)
  
  # Generate censoring times (e.g., dropout not due to events)
  participants <- generate_censoring_times(participants, config, scenario)
  
  # Determine observed time and status
  participants <- finalize_participant_data(participants, config)
  
  return(participants)
}

#' @title Generate enrollment times for participants
#' 
#' @param n_subjects Number of subjects
#' @param enrollment_rate Rate of enrollment per time unit
#' @return Vector of enrollment times
generate_enrollment_times <- function(n_subjects, enrollment_rate) {
  # Generate times based on exponential inter-arrival times
  if (enrollment_rate <= 0) {
    stop("Enrollment rate must be positive")
  }
  
  # Use exponential distribution for time between enrollments
  inter_arrival_times <- rexp(n_subjects, rate = enrollment_rate)
  
  # Cumulative sum to get actual enrollment times
  enrollment_times <- cumsum(inter_arrival_times)
  
  # Sort by enrollment time (first to last)
  enrollment_times <- sort(enrollment_times)
  
  return(enrollment_times)
}

#' @title Generate event times for participants based on arm-specific distributions
#' 
#' @param participants Participant data frame
#' @param config Trial configuration
#' @param scenario Simulation scenario
#' @return Updated participant data frame with event times
generate_event_times <- function(participants, config, scenario) {
  for (arm_name in names(config$arms)) {
    # Get arm-specific parameters
    arm_config <- config$arms[[arm_name]]
    arm_indices <- which(participants$arm == arm_name)
    n_arm <- length(arm_indices)
    
    # Modify parameters based on scenario
    tte_params <- adjust_params_for_scenario(arm_config$tte_params, arm_name, scenario)
    
    # Generate time-to-event based on specified distribution
    event_times <- switch(
      arm_config$tte_distribution,
      "exponential" = rexp(n_arm, rate = 1/tte_params$scale),
      "weibull" = rweibull(n_arm, shape = tte_params$shape, scale = tte_params$scale),
      "lognormal" = rlnorm(n_arm, meanlog = tte_params$meanlog, sdlog = tte_params$sdlog),
      "gamma" = rgamma(n_arm, shape = tte_params$shape, scale = tte_params$scale),
      # Default to Weibull if distribution is not recognized
      rweibull(n_arm, shape = 1.2, scale = 12)
    )
    
    # Assign to participant data frame
    participants$event_time[arm_indices] <- event_times
  }
  
  return(participants)
}

#' @title Adjust parameters based on simulation scenario
#' 
#' @param params Original parameters
#' @param arm_name Name of the arm
#' @param scenario Simulation scenario
#' @return Modified parameters for the scenario
adjust_params_for_scenario <- function(params, arm_name, scenario) {
  # For baseline scenario, return parameters unchanged
  if (scenario == "baseline") {
    return(params)
  }
  
  # For efficacy scenario, improve treatment arm parameters
  if (scenario == "efficacy" && arm_name == "treatment") {
    # For Weibull, increase scale to represent better efficacy
    if ("scale" %in% names(params)) {
      params$scale <- params$scale * 1.5  # Increase time to event by 50%
    }
    
    # For log-normal, decrease meanlog to represent better efficacy
    if ("meanlog" %in% names(params)) {
      params$meanlog <- params$meanlog + log(1.5)  # Increase time to event by 50%
    }
  }
  
  # For update scenario, parameters would typically be updated via Bayesian methods
  # This would be implemented in the apply_bayesian_updates function
  
  return(params)
}

#' @title Generate censoring times (dropouts not due to events)
#' 
#' @param participants Participant data frame
#' @param config Trial configuration
#' @param scenario Simulation scenario
#' @return Updated participant data frame with censoring times
generate_censoring_times <- function(participants, config, scenario) {
  n_subjects <- nrow(participants)
  
  # Base dropout rate from config
  dropout_rate <- config$dropout_rate
  
  # Adjust dropout rate based on scenario
  if (scenario == "efficacy") {
    # Lower dropout in efficacy scenario (better tolerability)
    dropout_rate <- dropout_rate * 0.8
  } else if (scenario == "update") {
    # Keep as is for update scenario or apply Bayesian update
  }
  
  # Generate censoring times using exponential distribution
  censoring_times <- rexp(n_subjects, rate = dropout_rate)
  
  # Add to enrollment time to get actual censoring time
  participants$censoring_time <- participants$enrollment_time + censoring_times
  
  return(participants)
}

#' @title Finalize participant data by determining observed time and event status
#' 
#' @param participants Participant data frame
#' @param config Trial configuration
#' @return Finalized participant data frame
finalize_participant_data <- function(participants, config) {
  # Calculate maximum possible follow-up time
  max_follow_up <- config$follow_up_duration
  
  # For each participant, determine the observed time and status
  for (i in 1:nrow(participants)) {
    enrollment_time <- participants$enrollment_time[i]
    event_time_raw <- participants$event_time[i]
    censoring_time_raw <- participants$censoring_time[i]
    
    # Calculate times relative to study start
    event_absolute <- enrollment_time + event_time_raw
    censoring_absolute <- enrollment_time + censoring_time_raw
    
    # End of follow-up
    end_of_follow_up <- enrollment_time + max_follow_up
    
    # Determine the earliest stopping time
    stopping_time <- min(event_absolute, censoring_absolute, end_of_follow_up)
    
    # Determine event status (1 = event, 0 = censored)
    event_status <- ifelse(stopping_time == event_absolute, 1, 0)
    
    # Determine retention status
    retention_status <- ifelse(stopping_time == censoring_absolute, "Dropout", 
                               ifelse(stopping_time == end_of_follow_up, "Completed", "Event"))
    
    # Calculate observed time (time from enrollment to stopping)
    observed_time <- stopping_time - enrollment_time
    
    # Update the participant data
    participants$observed_time[i] <- observed_time
    participants$event_status[i] <- event_status
    participants$retention_status[i] <- retention_status
  }
  
  return(participants)
}

#' @title Generate stratified participant data with additional covariates
#' 
#' @param config Trial configuration object
#' @param strata List of stratification factors and their distributions
#' @param covariate_effects List of covariate effects on outcomes
#' @param scenario Scenario type
#' @return Data frame of stratified participant data
generate_stratified_participants <- function(
    config,
    strata = list(
      site = list(
        levels = c("Site1", "Site2", "Site3", "Site4"),
        probs = c(0.3, 0.3, 0.2, 0.2)
      ),
      risk_group = list(
        levels = c("High", "Medium", "Low"),
        probs = c(0.25, 0.5, 0.25)
      )
    ),
    covariate_effects = list(
      risk_group = list(
        "High" = 1.5,    # Higher hazard (worse outcomes)
        "Medium" = 1.0,  # Reference
        "Low" = 0.7      # Lower hazard (better outcomes)
      ),
      age = 1.02  # 2% increase in hazard per year for age
    ),
    scenario = c("baseline", "efficacy", "update")
) {
  scenario <- match.arg(scenario)
  
  # First generate basic participant data
  participants <- generate_participants(config, scenario)
  
  # Add stratification factors
  for (factor_name in names(strata)) {
    factor_config <- strata[[factor_name]]
    
    participants[[factor_name]] <- sample(
      factor_config$levels,
      nrow(participants),
      replace = TRUE,
      prob = factor_config$probs
    )
  }
  
  # Apply covariate effects to event times
  if (!is.null(covariate_effects)) {
    for (i in 1:nrow(participants)) {
      hazard_multiplier <- 1.0
      
      # Apply categorical covariate effects
      for (covar_name in names(covariate_effects)) {
        if (covar_name %in% names(strata)) {
          # Categorical covariate
          level <- participants[[covar_name]][i]
          if (level %in% names(covariate_effects[[covar_name]])) {
            hazard_multiplier <- hazard_multiplier * covariate_effects[[covar_name]][[level]]
          }
        } else if (covar_name %in% names(participants)) {
          # Continuous covariate (like age)
          covar_value <- participants[[covar_name]][i]
          # Apply as multiplicative effect, e.g., for age: effect^(age-reference_age)
          reference_value <- 55  # Example reference age
          hazard_multiplier <- hazard_multiplier * 
            (covariate_effects[[covar_name]] ^ (covar_value - reference_value))
        }
      }
      
      # Adjust event time (higher hazard = shorter time)
      if (hazard_multiplier != 1.0) {
        participants$event_time[i] <- participants$event_time[i] / hazard_multiplier
      }
    }
    
    # Recalculate observed time and status after applying covariate effects
    participants <- finalize_participant_data(participants, config)
  }
  
  return(participants)
}

#' @title Generate a sample of synthetic patient data for testing
#' 
#' @param n_samples Number of patients to generate
#' @return A data frame with synthetic patient data
generate_test_data <- function(n_samples = 100) {
  # Create a simple test configuration
  test_config <- list(
    trial_name = "Test_Trial",
    n_subjects = n_samples,
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
    dropout_rate = 0.05,
    seed = 123
  )
  
  class(test_config) <- "trial_config"
  
  # Generate participants
  participants <- generate_participants(test_config, "baseline")
  
  return(participants)
}