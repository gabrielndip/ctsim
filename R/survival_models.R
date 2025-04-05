# Clinical Trial Simulation
# survival_models.R
# Created: 2025-03-31
# Description: Parametric and non-parametric survival models for clinical trial simulation

library(tidyverse)
library(survival)
library(flexsurv)
library(msm)

#' @title Fit Weibull survival model to trial data
#' 
#' @param data Participant data frame
#' @param formula Model formula 
#' @param strata Stratification variable
#' @return Fitted model and parameters
fit_weibull_model <- function(
    data,
    formula = Surv(observed_time, event_status) ~ arm,
    strata = NULL
) {
  # Convert to correct formula if strata is provided
  if (!is.null(strata)) {
    formula_str <- paste0("Surv(observed_time, event_status) ~ arm + strata(", strata, ")")
    formula <- as.formula(formula_str)
  }
  
  # Fit Weibull model using flexsurv
  model <- flexsurvreg(formula, data = data, dist = "weibull")
  
  # Extract parameters
  shape <- model$res["shape", "est"]
  scale_control <- exp(model$res["scale", "est"])
  
  # Calculate scale for treatment (if arm is in the model)
  if ("arm" %in% names(model$coefficients)) {
    scale_treatment <- exp(model$res["scale", "est"] + model$res["armtreatment", "est"])
  } else {
    scale_treatment <- scale_control
  }
  
  # Create parameter list
  params <- list(
    shape = shape,
    scale_control = scale_control,
    scale_treatment = scale_treatment,
    log_hr = -model$res["armtreatment", "est"] / shape,  # Log hazard ratio
    hazard_ratio = exp(-model$res["armtreatment", "est"] / shape)  # Hazard ratio
  )
  
  return(list(
    model = model,
    params = params
  ))
}

#' @title Fit exponential survival model to trial data
#' 
#' @param data Participant data frame
#' @param formula Model formula
#' @param strata Stratification variable
#' @return Fitted model and parameters
fit_exponential_model <- function(
    data,
    formula = Surv(observed_time, event_status) ~ arm,
    strata = NULL
) {
  # Convert to correct formula if strata is provided
  if (!is.null(strata)) {
    formula_str <- paste0("Surv(observed_time, event_status) ~ arm + strata(", strata, ")")
    formula <- as.formula(formula_str)
  }
  
  # Fit exponential model using flexsurv
  model <- flexsurvreg(formula, data = data, dist = "exponential")
  
  # Extract parameters (exponential is a special case of Weibull with shape=1)
  rate_control <- exp(-model$res["rate", "est"])
  
  # Calculate rate for treatment (if arm is in the model)
  if ("arm" %in% names(model$coefficients)) {
    rate_treatment <- exp(-model$res["rate", "est"] - model$res["armtreatment", "est"])
  } else {
    rate_treatment <- rate_control
  }
  
  # Create parameter list
  params <- list(
    rate_control = rate_control,
    rate_treatment = rate_treatment,
    log_hr = model$res["armtreatment", "est"],  # Log hazard ratio
    hazard_ratio = exp(model$res["armtreatment", "est"])  # Hazard ratio
  )
  
  return(list(
    model = model,
    params = params
  ))
}

#' @title Fit log-normal survival model to trial data
#' 
#' @param data Participant data frame
#' @param formula Model formula
#' @param strata Stratification variable
#' @return Fitted model and parameters
fit_lognormal_model <- function(
    data,
    formula = Surv(observed_time, event_status) ~ arm,
    strata = NULL
) {
  # Convert to correct formula if strata is provided
  if (!is.null(strata)) {
    formula_str <- paste0("Surv(observed_time, event_status) ~ arm + strata(", strata, ")")
    formula <- as.formula(formula_str)
  }
  
  # Fit log-normal model using flexsurv
  model <- flexsurvreg(formula, data = data, dist = "lnorm")
  
  # Extract parameters
  meanlog_control <- model$res["meanlog", "est"]
  sdlog <- model$res["sdlog", "est"]
  
  # Calculate meanlog for treatment (if arm is in the model)
  if ("arm" %in% names(model$coefficients)) {
    meanlog_treatment <- model$res["meanlog", "est"] + model$res["armtreatment", "est"]
    
    # Note: Hazard ratio in log-normal is not constant over time
    # We'll estimate it at the median survival time
    median_time <- exp(meanlog_control)
    hr_at_median <- dlnorm(median_time, meanlog_treatment, sdlog) / 
                   dlnorm(median_time, meanlog_control, sdlog) * 
                   (1 - plnorm(median_time, meanlog_control, sdlog)) / 
                   (1 - plnorm(median_time, meanlog_treatment, sdlog))
  } else {
    meanlog_treatment <- meanlog_control
    hr_at_median <- 1
  }
  
  # Create parameter list
  params <- list(
    meanlog_control = meanlog_control,
    meanlog_treatment = meanlog_treatment,
    sdlog = sdlog,
    hr_at_median = hr_at_median  # Approximate hazard ratio at median time
  )
  
  return(list(
    model = model,
    params = params
  ))
}

#' @title Fit Cox proportional hazards model
#' 
#' @param data Participant data frame
#' @param formula Model formula
#' @param strata Stratification variable
#' @return Fitted model and parameters
fit_cox_model <- function(
    data,
    formula = Surv(observed_time, event_status) ~ arm,
    strata = NULL
) {
  # Convert to correct formula if strata is provided
  if (!is.null(strata)) {
    formula_str <- paste0("Surv(observed_time, event_status) ~ arm + strata(", strata, ")")
    formula <- as.formula(formula_str)
  }
  
  # Fit Cox model
  model <- coxph(formula, data = data)
  
  # Extract hazard ratio and confidence interval
  if ("arm" %in% names(coef(model))) {
    hr <- exp(coef(model)["armtreatment"])
    hr_ci <- exp(confint(model)["armtreatment", ])
    p_value <- summary(model)$coefficients["armtreatment", "Pr(>|z|)"]
  } else {
    hr <- 1
    hr_ci <- c(NA, NA)
    p_value <- NA
  }
  
  # Create parameter list
  params <- list(
    hazard_ratio = hr,
    hr_ci_lower = hr_ci[1],
    hr_ci_upper = hr_ci[2],
    p_value = p_value
  )
  
  return(list(
    model = model,
    params = params
  ))
}

#' @title Compare different survival models
#' 
#' @param data Participant data frame
#' @param formula Base formula for models
#' @param time_points Time points for survival predictions
#' @return Comparison of models and AIC/BIC
compare_survival_models <- function(
    data,
    formula = Surv(observed_time, event_status) ~ arm,
    time_points = seq(0, max(data$observed_time), length.out = 20)
) {
  # Fit different models
  weibull_fit <- fit_weibull_model(data, formula)
  exponential_fit <- fit_exponential_model(data, formula)
  lognormal_fit <- fit_lognormal_model(data, formula)
  cox_fit <- fit_cox_model(data, formula)
  
  # Calculate AIC and BIC for parametric models
  model_metrics <- data.frame(
    Model = c("Weibull", "Exponential", "Log-normal"),
    AIC = c(AIC(weibull_fit$model), AIC(exponential_fit$model), AIC(lognormal_fit$model)),
    BIC = c(BIC(weibull_fit$model), BIC(exponential_fit$model), BIC(lognormal_fit$model))
  )
  
  # Sort by AIC
  model_metrics <- model_metrics[order(model_metrics$AIC), ]
  
  # Generate predictions for each model
  predictions <- data.frame(
    time = rep(time_points, 4),
    model = rep(c("Weibull", "Exponential", "Log-normal", "Cox"), each = length(time_points)),
    stringsAsFactors = FALSE
  )
  
  # Create newdata with treatment and control arms
  newdata_control <- data.frame(arm = "control")
  newdata_treatment <- data.frame(arm = "treatment")
  
  # Get predictions for control arm
  weibull_surv_control <- summary(weibull_fit$model, newdata = newdata_control, t = time_points)$survival
  exp_surv_control <- summary(exponential_fit$model, newdata = newdata_control, t = time_points)$survival
  lnorm_surv_control <- summary(lognormal_fit$model, newdata = newdata_control, t = time_points)$survival
  
  # Get predictions for treatment arm
  weibull_surv_treatment <- summary(weibull_fit$model, newdata = newdata_treatment, t = time_points)$survival
  exp_surv_treatment <- summary(exponential_fit$model, newdata = newdata_treatment, t = time_points)$survival
  lnorm_surv_treatment <- summary(lognormal_fit$model, newdata = newdata_treatment, t = time_points)$survival
  
  # Get Cox model predictions
  cox_surv <- survfit(cox_fit$model, newdata = data.frame(arm = c("control", "treatment")))
  cox_surv_times <- cox_surv$time
  cox_surv_control <- cox_surv$surv[, 1]
  cox_surv_treatment <- cox_surv$surv[, 2]
  
  # Interpolate Cox model predictions to match time points
  cox_interp_control <- approx(cox_surv_times, cox_surv_control, xout = time_points, rule = 2)$y
  cox_interp_treatment <- approx(cox_surv_times, cox_surv_treatment, xout = time_points, rule = 2)$y
  
  # Combine all predictions
  predictions$control <- c(weibull_surv_control, exp_surv_control, lnorm_surv_control, cox_interp_control)
  predictions$treatment <- c(weibull_surv_treatment, exp_surv_treatment, lnorm_surv_treatment, cox_interp_treatment)
  
  # Return model comparisons
  return(list(
    models = list(
      weibull = weibull_fit,
      exponential = exponential_fit,
      lognormal = lognormal_fit,
      cox = cox_fit
    ),
    metrics = model_metrics,
    predictions = predictions
  ))
}

#' @title Generate survival times from a specified parametric distribution
#' 
#' @param n Number of times to generate
#' @param distribution Distribution type
#' @param params Parameters for the distribution
#' @return Vector of survival times
generate_survival_times <- function(
    n,
    distribution = c("weibull", "exponential", "lognormal", "gamma"),
    params = list()
) {
  distribution <- match.arg(distribution)
  
  # Generate times based on specified distribution
  times <- switch(distribution,
    "weibull" = {
      if (!all(c("shape", "scale") %in% names(params))) {
        stop("Weibull distribution requires shape and scale parameters")
      }
      rweibull(n, shape = params$shape, scale = params$scale)
    },
    "exponential" = {
      if (!"rate" %in% names(params)) {
        stop("Exponential distribution requires rate parameter")
      }
      rexp(n, rate = params$rate)
    },
    "lognormal" = {
      if (!all(c("meanlog", "sdlog") %in% names(params))) {
        stop("Log-normal distribution requires meanlog and sdlog parameters")
      }
      rlnorm(n, meanlog = params$meanlog, sdlog = params$sdlog)
    },
    "gamma" = {
      if (!all(c("shape", "rate") %in% names(params))) {
        stop("Gamma distribution requires shape and rate parameters")
      }
      rgamma(n, shape = params$shape, rate = params$rate)
    }
  )
  
  return(times)
}

#' @title Create a competing risks survival model
#' 
#' @param cause_distributions List of distributions for each cause
#' @param cause_params List of parameters for each cause's distribution
#' @return Function that generates competing risks survival data
create_competing_risks_model <- function(
    cause_distributions = list(
      "primary" = "weibull",
      "secondary" = "exponential"
    ),
    cause_params = list(
      "primary" = list(shape = 1.2, scale = 12),
      "secondary" = list(rate = 1/24)
    )
) {
  # Validate inputs
  if (length(cause_distributions) != length(cause_params)) {
    stop("cause_distributions and cause_params must have the same length")
  }
  
  # Create the generator function
  generate_cr_data <- function(n) {
    # Initialize results
    times <- rep(Inf, n)
    causes <- rep(NA, n)
    
    # Generate times for each cause
    for (cause in names(cause_distributions)) {
      dist <- cause_distributions[[cause]]
      params <- cause_params[[cause]]
      
      # Generate times for this cause
      cause_times <- generate_survival_times(n, dist, params)
      
      # For each subject, if this cause happens earlier, update
      for (i in 1:n) {
        if (cause_times[i] < times[i]) {
          times[i] <- cause_times[i]
          causes[i] <- cause
        }
      }
    }
    
    return(data.frame(
      time = times,
      cause = causes
    ))
  }
  
  return(generate_cr_data)
}

#' @title Create a multi-state model for disease progression
#' 
#' @param states Vector of state names
#' @param transitions Matrix of transition intensities
#' @param arm_effect List of treatment effects on transitions
#' @return Function that generates multi-state model trajectories
create_multistate_model <- function(
    states = c("Healthy", "Disease", "Death"),
    transitions = matrix(c(
      0, 0.1, 0.05,
      0, 0, 0.2,
      0, 0, 0
    ), nrow = 3, byrow = TRUE),
    arm_effect = list(
      "control" = 1,
      "treatment" = 0.7  # 30% reduction in transition intensities
    )
) {
  # Create Q matrix
  qmatrix <- transitions
  dimnames(qmatrix) <- list(states, states)
  
  # Ensure diagonal elements are correct
  for (i in 1:nrow(qmatrix)) {
    qmatrix[i, i] <- -sum(qmatrix[i, -i])
  }
  
  # Create the generator function
  generate_msm_data <- function(n, arm, max_time = 24) {
    # Adjust transition intensities based on arm
    effect <- arm_effect[[arm]]
    if (is.null(effect)) effect <- 1
    
    q_adjusted <- qmatrix * effect
    
    # Generate trajectories
    data <- data.frame(
      id = 1:n,
      arm = arm,
      start_state = states[1],
      stringsAsFactors = FALSE
    )
    
    # For each subject, simulate trajectory
    trajectories <- list()
    
    for (i in 1:n) {
      # Start in first state
      current_state <- 1
      current_time <- 0
      trajectory <- data.frame(
        id = i,
        time = current_time,
        state = states[current_state],
        stringsAsFactors = FALSE
      )
      
      # Simulate until max_time or absorbing state
      while (current_time < max_time && sum(q_adjusted[current_state, ]) < 0) {
        # Calculate total transition intensity out of current state
        total_intensity <- -q_adjusted[current_state, current_state]
        
        # If zero, we're in an absorbing state
        if (total_intensity == 0) break
        
        # Generate time to next transition
        time_to_next <- rexp(1, rate = total_intensity)
        current_time <- current_time + time_to_next
        
        # If beyond max_time, stop
        if (current_time > max_time) break
        
        # Determine next state
        transition_probs <- q_adjusted[current_state, ] / total_intensity
        transition_probs[current_state] <- 0  # No self-transitions
        next_state <- sample(1:length(states), 1, prob = transition_probs)
        
        # Add to trajectory
        trajectory <- rbind(trajectory, data.frame(
          id = i,
          time = current_time,
          state = states[next_state],
          stringsAsFactors = FALSE
        ))
        
        current_state <- next_state
      }
      
      trajectories[[i]] <- trajectory
    }
    
    # Combine all trajectories
    result <- do.call(rbind, trajectories)
    return(result)
  }
  
  return(generate_msm_data)
}

