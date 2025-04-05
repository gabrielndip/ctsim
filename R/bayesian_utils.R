# Clinical Trial Simulation
# bayesian_utils.R
# Created: 2025-03-31
# Description: Bayesian utilities for clinical trial simulation

library(tidyverse)
library(brms)
library(rstanarm)

#' @title Create log hazard ratio prior distribution
#' 
#' @param distribution Prior distribution type
#' @param mean Mean of the prior distribution
#' @param sd Standard deviation of the prior distribution
#' @param alpha Shape parameter for beta or gamma distributions
#' @param beta Scale parameter for beta or gamma distributions
#' @return A list containing prior distribution specifications
create_hr_prior <- function(
    distribution = c("normal", "student_t", "cauchy", "beta", "gamma"),
    mean = 0,
    sd = 1,
    alpha = NULL,
    beta = NULL
) {
  distribution <- match.arg(distribution)
  
  prior <- list(
    distribution = distribution,
    parameters = switch(distribution,
                        normal = list(mean = mean, sd = sd),
                        student_t = list(df = 3, mean = mean, sd = sd),
                        cauchy = list(location = mean, scale = sd),
                        beta = list(alpha = alpha, beta = beta),
                        gamma = list(shape = alpha, rate = beta)),
    is_informative = !(distribution == "normal" && mean == 0 && sd > 10)
  )
  
  class(prior) <- "bayesian_prior"
  return(prior)
}

#' @title Update Bayesian prior with new data
#' 
#' @param prior Bayesian prior object
#' @param data New data for updating
#' @param formula Formula for the model
#' @param additional_args Additional arguments for the model
#' @return Updated posterior distribution
update_prior <- function(
    prior,
    data,
    formula = Surv(observed_time, event_status) ~ arm,
    additional_args = list()
) {
  # Set up prior specs for brms
  if (prior$distribution == "normal") {
    prior_spec <- prior(normal(prior$parameters$mean, prior$parameters$sd), class = "b", coef = "armtreatment")
  } else if (prior$distribution == "student_t") {
    prior_spec <- prior(student_t(prior$parameters$df, prior$parameters$mean, prior$parameters$sd), 
                         class = "b", coef = "armtreatment")
  } else if (prior$distribution == "cauchy") {
    prior_spec <- prior(cauchy(prior$parameters$location, prior$parameters$scale), 
                         class = "b", coef = "armtreatment")
  } else {
    # For beta and gamma, we may need to transform or use a different approach
    # For simplicity, we'll use a normal approximation
    prior_spec <- prior(normal(0, 10), class = "b", coef = "armtreatment")
  }
  
  # Fit Bayesian model using brms
  fit_args <- c(
    list(
      formula = formula,
      data = data,
      family = brmsfamily("cox"),
      prior = prior_spec,
      chains = 4,
      iter = 2000,
      warmup = 1000,
      cores = 2,
      seed = 123
    ),
    additional_args
  )
  
  model <- do.call(brm, fit_args)
  
  # Extract posterior distribution
  posterior_samples <- posterior_samples(model)
  
  # Calculate posterior summary
  posterior <- list(
    model = model,
    samples = posterior_samples,
    mean = mean(posterior_samples$b_armtreatment),
    median = median(posterior_samples$b_armtreatment),
    sd = sd(posterior_samples$b_armtreatment),
    quantiles = quantile(posterior_samples$b_armtreatment, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)),
    prior = prior
  )
  
  class(posterior) <- "bayesian_posterior"
  return(posterior)
}

#' @title Calculate probability of treatment effect exceeding threshold
#' 
#' @param posterior Bayesian posterior object
#' @param threshold Effect size threshold (log hazard ratio scale)
#' @return Probability that effect exceeds threshold
calculate_effect_probability <- function(posterior, threshold = log(0.8)) {
  # For Cox models, negative values indicate reduced hazard (benefit)
  # So we want to calculate P(log_hr < threshold)
  prob <- mean(posterior$samples$b_armtreatment < threshold)
  return(prob)
}

#' @title Generate predictions from Bayesian survival model
#' 
#' @param posterior Bayesian posterior object
#' @param newdata Data to generate predictions for
#' @param times Time points to predict at
#' @return Predicted survival probabilities
predict_survival <- function(posterior, newdata, times = seq(0, 24, by = 0.5)) {
  # Generate survival predictions from the model
  predictions <- predict(posterior$model, newdata = newdata, type = "response")
  
  # Format predictions for return
  pred_df <- data.frame(
    time = rep(times, each = nrow(newdata)),
    arm = rep(newdata$arm, length(times)),
    survival = as.vector(predictions)
  )
  
  return(pred_df)
}

#' @title Implement adaptive randomization based on posterior probabilities
#' 
#' @param posterior Bayesian posterior object
#' @param n_subjects Number of subjects to randomize
#' @param min_allocation Minimum allocation proportion for any arm
#' @return Vector of arm assignments
adaptive_randomization <- function(
    posterior,
    n_subjects = 100,
    min_allocation = 0.2
) {
  # Calculate probability that treatment is better than control
  prob_better <- calculate_effect_probability(posterior)
  
  # Calculate allocation ratio for treatment
  # More evidence of benefit = higher allocation to treatment
  treatment_alloc <- min(max(prob_better, min_allocation), 1 - min_allocation)
  control_alloc <- 1 - treatment_alloc
  
  # Generate random assignments
  arms <- sample(c("control", "treatment"), n_subjects, replace = TRUE, 
                 prob = c(control_alloc, treatment_alloc))
  
  return(list(
    arms = arms,
    allocation = c(control = control_alloc, treatment = treatment_alloc),
    prob_better = prob_better
  ))
}

