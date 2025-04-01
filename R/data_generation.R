# Clinical Trial Simulation
# data_generation.R
# Created: 2025-04-01
# Description: Functions for generating parameter configurations for simulated clinical trial data

#' @title Create a trial parameters data structure 
#' @param trial_name Character string for trial identifier
#' @param n_subjects Number of subjects to enroll
#' @param enrollment_rate Number of subjects enrolled per time unit
#' @param follow_up_duration Duration of follow-up in time units
#' @param arms List of trial arms with their parameters
#' @param events_required Number of events required for analysis
#' @param dropout_rate Baseline dropout rate (non-event related)
#' @param seed Random seed for reproducibility
#' @param prior_data Optional list containing prior data for Bayesian updates
#' @return A list containing the trial configuration
create_trial_config <- function(
    trial_name = "Trial_01",
    n_subjects = 100,
    enrollment_rate = 2,  # subjects per month as time unit
    follow_up_duration = 24,  # in months as time unit
    arms = list(
      control = list(
        name = "Control",
        allocation = 0.5,  # proportion of subjects
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
    dropout_rate = 0.1,  # per month as time unit
    seed = 2202,
    prior_data = NULL
) {
  # Validate inputs???
  if (sum(sapply(arms, function(x) x$allocation)) != 1) {
    stop("Arm allocations must sum to 1")
  }
  
  # Configure the trial
  config <- list(
    trial_name = trial_name,
    n_subjects = n_subjects,
    enrollment_rate = enrollment_rate,
    follow_up_duration = follow_up_duration,
    arms = arms,
    events_required = events_required,
    dropout_rate = dropout_rate,
    seed = seed,
    prior_data = prior_data,
    simulation_time = NULL,  # Will be calculated during simulation
    interim_analyses = list()
  )
  
  class(config) <- "trial_config"
  return(config)
}

#' @title Add interim analysis to trial configuration
#' 
#' @param config Trial configuration object
#' @param timing When to conduct analysis (proportion of events or time)
#' @param type Type of interim analysis ("events" or "time")
#' @param decision_rules List of rules for trial decisions
#' @return Updated trial configuration object
add_interim_analysis <- function(
    config,
    timing,
    type = c("events", "time"),
    decision_rules = list(
      efficacy = function(results) { results$hazard_ratio < 0.7 & results$p_value < 0.05 },
      futility = function(results) { results$hazard_ratio > 0.9 | results$p_value > 0.3 }
    )
) {
  type <- match.arg(type)
  
  interim <- list(
    timing = timing,
    type = type,
    decision_rules = decision_rules
  )
  
  config$interim_analyses <- c(config$interim_analyses, list(interim))
  return(config)
}

#' @title Structure for holding simulation results
#' 
#' @param config Trial configuration used
#' @return A structured list for holding simulation results
create_simulation_results <- function(config) {
  results <- list(
    config = config,
    participants = data.frame(),
    events = data.frame(),
    analyses = list(),
    trial_completion_time = NA,
    event_count = 0,
    final_analysis = list(),
    scenario = "baseline",
    bayesian_updates = list()
  )
  
  class(results) <- "trial_simulation"
  return(results)
}

#' @trial Create a Bayesian prior specification
#' 
#' @param distribution Prior distribution type (e.g., "normal", "gamma")
#' @param parameters List of parameters for the distribution
#' @return A structured prior specification
create_bayesian_prior <- function(
    distribution = c("normal", "gamma", "beta", "informative", "non_informative"),
    parameters = list()
) {
  distribution <- match.arg(distribution)
  
  prior <- list(
    distribution = distribution,
    parameters = parameters
  )
  
  class(prior) <- "bayesian_prior"
  return(prior)
}

#' @example Create a trial configuration with Bayesian priors
example_trial_with_priors <- function() {
  # Create treatment effect prior
  hr_prior <- create_bayesian_prior(
    distribution = "normal",
    parameters = list(mean = 0.8, sd = 0.2)
  )
  
  # Create dropout rate prior
  dropout_prior <- create_bayesian_prior(
    distribution = "beta",
    parameters = list(shape1 = 2, shape2 = 20)
  )
  
  # Create trial configuration
  config <- create_trial_config(
    trial_name = "BayesianTrial_01",
    n_subjects = 200,
    arms = list(
      control = list(
        name = "Control",
        allocation = 0.5,
        tte_distribution = "weibull",
        tte_params = list(shape = 1.2, scale = 12),
        prior = NULL
      ),
      treatment = list(
        name = "Treatment",
        allocation = 0.5,
        tte_distribution = "weibull",
        tte_params = list(shape = 1.2, scale = 18),
        prior = hr_prior
      )
    ),
    dropout_rate = 0.1,
    prior_data = list(
      dropout_prior = dropout_prior,
      previous_trials = NULL  # Could include data from previous trials
    )
  )
  
  return(config)
}

