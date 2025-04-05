# Common survival analysis functions for ctsim package
# This file contains helper functions used across the package

#' Create a Kaplan-Meier plot from survival data
#'
#' @param time Vector of observed times
#' @param event Vector of event indicators (1=event, 0=censored)
#' @param group Vector of group assignments
#' @param data Data frame containing the variables
#' @return A ggplot object with the Kaplan-Meier plot
km_plot <- function(time, event, group, data) {
  # Fit Kaplan-Meier curve
  require(survival)
  require(ggplot2)
  
  km_fit <- survfit(Surv(time, event) ~ group, data=data)
  
  # Create plot data frame
  plot_data <- data.frame(
    time = km_fit$time,
    surv = km_fit$surv,
    strata = factor(km_fit$strata)
  )
  
  # Create plot
  ggplot(plot_data, aes(x = time, y = surv, color = strata)) +
    geom_step() +
    ylim(0, 1) +
    labs(x = "Time", y = "Survival Probability", 
         title = "Kaplan-Meier Survival Estimate") +
    theme_minimal()
}

#' Fit parametric survival models
#'
#' @param time Vector of observed times
#' @param event Vector of event indicators (1=event, 0=censored)
#' @param covariates Matrix or data frame of covariates
#' @param data Data frame containing all variables
#' @return List of fitted parametric models
fit_parametric_models <- function(time, event, covariates, data) {
  require(survival)
  require(flexsurv)
  
  # Create formula
  formula <- as.formula(paste("Surv(time, event) ~", paste(covariates, collapse = " + ")))
  
  # Fit models
  weib_model <- survreg(formula, data=data, dist="weibull")
  exp_model <- survreg(formula, data=data, dist="exponential")
  logn_model <- survreg(formula, data=data, dist="lognormal")
  
  # Return as list
  list(
    weibull = weib_model,
    exponential = exp_model,
    lognormal = logn_model
  )
}

#' Generate survival times from a Weibull model
#'
#' @param n Number of times to generate
#' @param covariates Matrix of covariate values
#' @param model Fitted survreg model for Weibull distribution
#' @return Vector of generated survival times
generate_weibull_times <- function(n, covariates, model) {
  # Calculate linear predictor
  linear_pred <- model$coefficients[1] + 
    as.matrix(covariates) %*% model$coefficients[-1]
  
  # Convert to Weibull parameters
  scale <- exp(linear_pred)
  shape <- 1/model$scale
  
  # Generate random survival times
  rweibull(n, shape=shape, scale=scale)
}

#' Create observed survival data with censoring
#'
#' @param event_times Vector of true event times
#' @param max_follow_up Maximum follow-up time
#' @param dropout_rate Rate parameter for dropout process
#' @return List with observed times and event indicators
create_observed_data <- function(event_times, max_follow_up, dropout_rate=0.1) {
  n <- length(event_times)
  
  # Generate random dropout times
  dropout_times <- rexp(n, rate=dropout_rate)
  
  # Administrative censoring at max follow-up
  censor_times <- pmin(dropout_times, max_follow_up)
  
  # Determine observed time (minimum of event time and censoring time)
  observed_times <- pmin(event_times, censor_times)
  
  # Create event indicators (1 if event observed, 0 if censored)
  event_indicators <- as.numeric(event_times <= censor_times)
  
  # Return as list
  list(
    time = observed_times,
    status = event_indicators
  )
}