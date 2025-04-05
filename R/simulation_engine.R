# Clinical Trial Simulation
# simulation_engine.R
# Created: 2025-04-01
# Description: Simulation Engine Core Functions
# This file contains the core functions for running clinical trial simulations

library(survival)
library(dplyr)
library(ggplot2)

#' @title Run a full clinical trial simulation
#' 
#' @param config Trial configuration object created by create_trial_config()
#' @param n_sims Number of simulation iterations
#' @param scenario Scenario to run ("baseline", "efficacy", "update")
#' @param verbose Logical, whether to print progress information
#' @return List containing simulation results
run_trial_simulation <- function(
    config,
    n_sims = 1,
    scenario = c("baseline", "efficacy", "update"),
    verbose = TRUE
) {
  scenario <- match.arg(scenario)
  
  if (verbose) {
    cat("Starting simulation of trial:", config$trial_name, "\n")
    cat("Scenario:", scenario, "\n")
  }
  
  # Set random seed for reproducibility
  set.seed(config$seed)
  
  # Initialize results container
  all_results <- list()
  
  # Run multiple simulations if requested
  for (sim in 1:n_sims) {
    if (verbose && n_sims > 1) {
      cat("Running simulation", sim, "of", n_sims, "\n")
    }
    
    # Initialize single simulation result
    sim_result <- create_simulation_results(config)
    sim_result$scenario <- scenario
    
    # Generate participant data
    sim_result$participants <- generate_participants(config, scenario)
    
    # Simulate enrollment and event times
    sim_result <- simulate_trial_events(sim_result)
    
    # Run interim analyses if specified
    sim_result <- run_interim_analyses(sim_result)
    
    # Run final analysis
    sim_result <- run_final_analysis(sim_result)
    
    # If using Bayesian updates, apply them
    if (scenario == "update" && !is.null(config$prior_data)) {
      sim_result <- apply_bayesian_updates(sim_result)
    }
    
    # Store results for this simulation
    all_results[[sim]] <- sim_result
  }
  
  # If only one simulation was run, return it directly
  if (n_sims == 1) {
    return(all_results[[1]])
  }
  
  # Otherwise, summarize results across simulations
  summary_results <- summarize_simulations(all_results)
  return(list(individual_sims = all_results, summary = summary_results))
}

#' @title Simulate trial events (enrollment, time-to-event, dropouts)
#' 
#' @param sim_result Simulation results object
#' @return Updated simulation results with events data
simulate_trial_events <- function(sim_result) {
  # Extract participant data
  participants <- sim_result$participants
  config <- sim_result$config
  
  # Create a sorted timeline of all events
  event_timeline <- data.frame(
    participant_id = character(),
    time = numeric(),
    event_type = character(),
    arm = character(),
    stringsAsFactors = FALSE
  )
  
  # Add enrollment events to timeline
  enrollment_events <- data.frame(
    participant_id = participants$participant_id,
    time = participants$enrollment_time,
    event_type = "enrollment",
    arm = participants$arm,
    stringsAsFactors = FALSE
  )
  event_timeline <- rbind(event_timeline, enrollment_events)
  
  # Add observed events (event or censoring) to timeline
  observed_events <- data.frame(
    participant_id = participants$participant_id,
    time = participants$enrollment_time + participants$observed_time,
    event_type = ifelse(participants$event_status == 1, "event", "censoring"),
    arm = participants$arm,
    stringsAsFactors = FALSE
  )
  event_timeline <- rbind(event_timeline, observed_events)
  
  # Sort by time
  event_timeline <- event_timeline[order(event_timeline$time), ]
  
  # Initialize counters for tracking trial progress
  enrolled_count <- 0
  event_count <- 0
  censored_count <- 0
  arm_counts <- setNames(rep(0, length(names(config$arms))), names(config$arms))
  arm_event_counts <- setNames(rep(0, length(names(config$arms))), names(config$arms))
  
  # Initialize detailed event log
  event_log <- data.frame(
    time = numeric(),
    participant_id = character(),
    event_type = character(),
    arm = character(),
    enrolled_count = numeric(),
    event_count = numeric(),
    censored_count = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Process events in chronological order
  for (i in 1:nrow(event_timeline)) {
    event <- event_timeline[i, ]
    
    # Update counters based on event type
    if (event$event_type == "enrollment") {
      enrolled_count <- enrolled_count + 1
      arm_counts[event$arm] <- arm_counts[event$arm] + 1
    } else if (event$event_type == "event") {
      event_count <- event_count + 1
      arm_event_counts[event$arm] <- arm_event_counts[event$arm] + 1
    } else if (event$event_type == "censoring") {
      censored_count <- censored_count + 1
    }
    
    # Add to event log
    event_log <- rbind(event_log, data.frame(
      time = event$time,
      participant_id = event$participant_id,
      event_type = event$event_type,
      arm = event$arm,
      enrolled_count = enrolled_count,
      event_count = event_count,
      censored_count = censored_count,
      stringsAsFactors = FALSE
    ))
    
    # Check if we've reached the required number of events for trial completion
    if (event_count >= config$events_required) {
      sim_result$trial_completion_time <- event$time
      break
    }
  }
  
  # If trial did not reach required events, set completion time to the last event time
  if (is.na(sim_result$trial_completion_time) && nrow(event_log) > 0) {
    sim_result$trial_completion_time <- max(event_log$time)
  }
  
  # Update simulation results
  sim_result$events <- event_log
  sim_result$event_count <- event_count
  
  # Calculate additional metrics
  sim_result$metrics <- list(
    enrollment_duration = ifelse(enrolled_count == config$n_subjects,
                               max(participants$enrollment_time), NA),
    study_duration = sim_result$trial_completion_time,
    event_rate = event_count / enrolled_count,
    censoring_rate = censored_count / enrolled_count,
    arm_enrollment = arm_counts,
    arm_events = arm_event_counts
  )
  
  return(sim_result)
}

#' @title Run interim analyses as specified in the config
#' 
#' @param sim_result Simulation results object
#' @return Updated simulation results with interim analyses results
run_interim_analyses <- function(sim_result) {
  # Extract configuration and event data
  config <- sim_result$config
  events <- sim_result$events
  participants <- sim_result$participants
  
  # If no interim analyses are defined, return unchanged
  if (length(config$interim_analyses) == 0) {
    sim_result$analyses <- list()
    return(sim_result)
  }
  
  # Initialize list to store interim analysis results
  interim_results <- list()
  
  # Process each interim analysis
  for (i in seq_along(config$interim_analyses)) {
    interim_config <- config$interim_analyses[[i]]
    
    # Determine when this interim analysis should occur
    if (interim_config$type == "events") {
      # Based on number of events
      event_threshold <- round(interim_config$timing * config$events_required)
      
      # Find the timepoint when we reach this number of events
      event_rows <- which(events$event_type == "event")
      if (length(event_rows) >= event_threshold) {
        analysis_time <- events$time[event_rows[event_threshold]]
        analysis_triggered <- TRUE
      } else {
        analysis_triggered <- FALSE
      }
    } else if (interim_config$type == "time") {
      # Based on calendar time
      analysis_time <- interim_config$timing
      analysis_triggered <- max(events$time) >= analysis_time
    }
    
    # Skip if analysis wasn't triggered
    if (!analysis_triggered) {
      next
    }
    
    # Select data up to the analysis timepoint
    analysis_events <- events[events$time <= analysis_time, ]
    current_participants <- participants
    
    # Update event status and observed time for participants based on interim analysis time
    for (p in 1:nrow(current_participants)) {
      enrollment_time <- current_participants$enrollment_time[p]
      event_time_absolute <- enrollment_time + current_participants$event_time[p]
      censoring_time_absolute <- enrollment_time + current_participants$censoring_time[p]
      
      # If the event or censoring would occur after analysis time, censor at analysis time
      if (event_time_absolute > analysis_time && censoring_time_absolute > analysis_time) {
        current_participants$event_status[p] <- 0  # Censored
        current_participants$observed_time[p] <- analysis_time - enrollment_time
        current_participants$retention_status[p] <- "Ongoing"
      }
    }
    
    # Run survival analysis on interim data
    analysis_result <- try({
      # Fit Cox proportional hazards model
      survival_formula <- Surv(observed_time, event_status) ~ arm
      cox_model <- coxph(survival_formula, data = current_participants)
      
      # Extract hazard ratio (treatment vs. control)
      hr <- exp(coef(cox_model)["armtreatment"])
      hr_ci <- exp(confint(cox_model)["armtreatment", ])
      p_value <- summary(cox_model)$coefficients["armtreatment", "Pr(>|z|)"]
      
      # Calculate event rates by arm
      arm_event_rates <- tapply(current_participants$event_status, 
                                current_participants$arm, 
                                function(x) sum(x) / length(x))
      
      list(
        time = analysis_time,
        event_count = sum(current_participants$event_status),
        enrolled_count = nrow(current_participants),
        cox_model = cox_model,
        hazard_ratio = hr,
        hr_ci_lower = hr_ci[1],
        hr_ci_upper = hr_ci[2],
        p_value = p_value,
        arm_event_rates = arm_event_rates
      )
    }, silent = TRUE)
    
    # Handle any errors in analysis
    if (inherits(analysis_result, "try-error")) {
      analysis_result <- list(
        time = analysis_time,
        event_count = sum(current_participants$event_status),
        enrolled_count = nrow(current_participants),
        error = TRUE,
        error_message = attr(analysis_result, "condition")$message
      )
    }
    
    # Apply decision rules
    if (!inherits(analysis_result, "try-error") && !is.null(interim_config$decision_rules)) {
      # Apply efficacy and futility rules
      if ("efficacy" %in% names(interim_config$decision_rules)) {
        analysis_result$stop_for_efficacy <- interim_config$decision_rules$efficacy(analysis_result)
      } else {
        analysis_result$stop_for_efficacy <- FALSE
      }
      
      if ("futility" %in% names(interim_config$decision_rules)) {
        analysis_result$stop_for_futility <- interim_config$decision_rules$futility(analysis_result)
      } else {
        analysis_result$stop_for_futility <- FALSE
      }
      
      # Determine if trial should stop
      analysis_result$stop_trial <- analysis_result$stop_for_efficacy || analysis_result$stop_for_futility
    } else {
      analysis_result$stop_for_efficacy <- FALSE
      analysis_result$stop_for_futility <- FALSE
      analysis_result$stop_trial <- FALSE
    }
    
    # Store the interim analysis result
    interim_results[[i]] <- analysis_result
    
    # If decision is to stop the trial, break out of the loop
    if (analysis_result$stop_trial) {
      sim_result$trial_stopped_early <- TRUE
      sim_result$trial_completion_time <- analysis_time
      sim_result$stop_reason <- ifelse(analysis_result$stop_for_efficacy, "efficacy", "futility")
      break
    }
  }
  
  # Update simulation results
  sim_result$analyses <- interim_results
  if (length(interim_results) > 0) {
    sim_result$interim_decisions <- sapply(interim_results, function(x) x$stop_trial)
  }
  
  return(sim_result)
}

#' @title Run the final analysis for the simulated trial
#' 
#' @param sim_result Simulation results object
#' @return Updated simulation results with final analysis
run_final_analysis <- function(sim_result) {
  # Check if trial was stopped early at an interim analysis
  if (!is.null(sim_result$trial_stopped_early) && sim_result$trial_stopped_early) {
    # Use the last interim analysis as the final analysis
    last_interim <- sim_result$analyses[[length(sim_result$analyses)]]
    sim_result$final_analysis <- last_interim
    sim_result$final_analysis$is_interim <- TRUE
    return(sim_result)
  }
  
  # Extract participant data at trial completion
  participants <- sim_result$participants
  config <- sim_result$config
  
  # If trial completed with the required number of events, run final analysis
  try({
    # Fit Cox proportional hazards model
    survival_formula <- Surv(observed_time, event_status) ~ arm
    cox_model <- coxph(survival_formula, data = participants)
    
    # Extract hazard ratio (treatment vs. control)
    hr <- exp(coef(cox_model)["armtreatment"])
    hr_ci <- exp(confint(cox_model)["armtreatment", ])
    p_value <- summary(cox_model)$coefficients["armtreatment", "Pr(>|z|)"]
    
    # Conduct log-rank test
    surv_diff <- survdiff(survival_formula, data = participants)
    log_rank_p <- 1 - pchisq(surv_diff$chisq, df = length(surv_diff$n) - 1)
    
    # Calculate Kaplan-Meier estimates for plotting
    km_fit <- survfit(Surv(observed_time, event_status) ~ arm, data = participants)
    
    # Calculate median survival time by arm
    median_survival <- summary(km_fit)$table[, "median"]
    
    # Calculate event rates by arm
    arm_event_counts <- tapply(participants$event_status, participants$arm, sum)
    arm_totals <- table(participants$arm)
    arm_event_rates <- arm_event_counts / arm_totals
    
    # Perform parametric survival analysis using Weibull model
    weibull_model <- try(survreg(survival_formula, data = participants, dist = "weibull"), silent = TRUE)
    if (!inherits(weibull_model, "try-error")) {
      # Convert to Weibull parameters in the usual parameterization
      weibull_shape <- 1 / weibull_model$scale
      weibull_scale_control <- exp(coef(weibull_model)[1])
      weibull_scale_treatment <- exp(coef(weibull_model)[1] + coef(weibull_model)[2])
      
      parametric_models <- list(
        weibull = list(
          model = weibull_model,
          shape = weibull_shape,
          scale_control = weibull_scale_control,
          scale_treatment = weibull_scale_treatment
        )
      )
    } else {
      parametric_models <- list(
        error = TRUE,
        error_message = attr(weibull_model, "condition")$message
      )
    }
    
    # Store final analysis results
    sim_result$final_analysis <- list(
      time = sim_result$trial_completion_time,
      event_count = sum(participants$event_status),
      enrolled_count = nrow(participants),
      cox_model = cox_model,
      hazard_ratio = hr,
      hr_ci_lower = hr_ci[1],
      hr_ci_upper = hr_ci[2],
      p_value = p_value,
      log_rank_p = log_rank_p,
      km_fit = km_fit,
      median_survival = median_survival,
      arm_event_rates = arm_event_rates,
      parametric_models = parametric_models,
      is_interim = FALSE
    )
    
    # Determine trial outcome (success/failure)
    sim_result$success <- p_value < 0.05  # Default threshold
  }, silent = TRUE)
  
  # Handle errors in analysis
  if (is.null(sim_result$final_analysis) || !is.list(sim_result$final_analysis)) {
    sim_result$final_analysis <- list(
      time = sim_result$trial_completion_time,
      event_count = sum(participants$event_status),
      enrolled_count = nrow(participants),
      error = TRUE,
      error_message = "Error in final analysis"
    )
    sim_result$success <- FALSE
  }
  
  return(sim_result)
}

#' @title Apply Bayesian updates based on new data
#' 
#' @param sim_result Simulation results object
#' @return Updated simulation results with Bayesian posterior updates
apply_bayesian_updates <- function(sim_result) {
  # Extract configuration and data
  config <- sim_result$config
  participants <- sim_result$participants
  
  # Check if prior data is available
  if (is.null(config$prior_data)) {
    sim_result$bayesian_updates <- list(
      message = "No prior data provided for Bayesian updates."
    )
    return(sim_result)
  }
  
  # Initialize Bayesian updates list
  bayesian_updates <- list()
  
  # Extract treatment arm prior if available
  if (!is.null(config$arms$treatment$prior) && 
      inherits(config$arms$treatment$prior, "bayesian_prior")) {
    
    # Get treatment effect prior
    hr_prior <- config$arms$treatment$prior
    
    try({
      # Update prior with current data
      hr_posterior <- update_prior(
        prior = hr_prior,
        data = participants,
        formula = Surv(observed_time, event_status) ~ arm
      )
      
      # Calculate probability of efficacy
      prob_efficacy <- calculate_effect_probability(hr_posterior, threshold = log(0.8))
      
      # Store posterior distribution
      bayesian_updates$hr_posterior <- hr_posterior
      bayesian_updates$prob_efficacy <- prob_efficacy
      
      # Generate predictive survival curves
      newdata <- data.frame(
        arm = c("control", "treatment"),
        stringsAsFactors = FALSE
      )
      predicted_survival <- predict_survival(hr_posterior, newdata)
      bayesian_updates$predicted_survival <- predicted_survival
      
      # Calculate probability of superiority
      bayesian_updates$prob_superiority <- calculate_effect_probability(hr_posterior, threshold = 0)
      
      # Simulate adaptive randomization for future subjects
      if (config$n_subjects > nrow(participants)) {
        future_n <- config$n_subjects - nrow(participants)
        adaptive_rand <- adaptive_randomization(hr_posterior, future_n)
        bayesian_updates$adaptive_randomization <- adaptive_rand
      }
    }, silent = TRUE)
  }
  
  # Check for dropout rate prior
  if (!is.null(config$prior_data$dropout_prior)) {
    # Placeholder for dropout rate update
    # This would be implemented with beta-binomial conjugate analysis
    bayesian_updates$dropout_update <- list(
      message = "Dropout rate update would be implemented here."
    )
  }
  
  # If Bayesian updates failed, return a message
  if (length(bayesian_updates) == 0) {
    bayesian_updates <- list(
      message = "Bayesian updates failed or no valid priors were found."
    )
  }
  
  # Update simulation results
  sim_result$bayesian_updates <- bayesian_updates
  
  return(sim_result)
}

#' @title Summarize results across multiple simulations
#' 
#' @param all_results List of simulation results
#' @param metrics Metrics to include in summary
#' @return Summary statistics across simulations
summarize_simulations <- function(
    all_results,
    metrics = c("success_rate", "trial_duration", "hr", "events", "early_stopping")
) {
  # Ensure all_results is a list
  if (!is.list(all_results)) {
    stop("all_results must be a list of simulation results")
  }
  
  # Number of simulations
  n_sims <- length(all_results)
  
  # Initialize summary statistics
  summary <- list(
    n_sims = n_sims,
    config = all_results[[1]]$config, # Use config from first simulation
    metrics = metrics
  )
  
  # Create vectors to store key metrics
  success <- numeric(n_sims)
  trial_duration <- numeric(n_sims)
  hazard_ratios <- numeric(n_sims)
  hr_ci_lower <- numeric(n_sims)
  hr_ci_upper <- numeric(n_sims)
  p_values <- numeric(n_sims)
  event_counts <- numeric(n_sims)
  early_stopping <- logical(n_sims)
  early_stopping_reason <- character(n_sims)
  
  # Extract metrics from each simulation
  for (i in 1:n_sims) {
    sim <- all_results[[i]]
    
    # Trial success
    success[i] <- ifelse(is.null(sim$success), NA, sim$success)
    
    # Trial duration
    trial_duration[i] <- ifelse(is.null(sim$trial_completion_time), NA, sim$trial_completion_time)
    
    # Hazard ratio and confidence interval
    if (!is.null(sim$final_analysis) && !is.null(sim$final_analysis$hazard_ratio)) {
      hazard_ratios[i] <- sim$final_analysis$hazard_ratio
      hr_ci_lower[i] <- sim$final_analysis$hr_ci_lower
      hr_ci_upper[i] <- sim$final_analysis$hr_ci_upper
      p_values[i] <- sim$final_analysis$p_value
    } else {
      hazard_ratios[i] <- NA
      hr_ci_lower[i] <- NA
      hr_ci_upper[i] <- NA
      p_values[i] <- NA
    }
    
    # Event counts
    event_counts[i] <- sim$event_count
    
    # Early stopping
    early_stopping[i] <- ifelse(is.null(sim$trial_stopped_early), FALSE, sim$trial_stopped_early)
    if (early_stopping[i] && !is.null(sim$stop_reason)) {
      early_stopping_reason[i] <- sim$stop_reason
    } else {
      early_stopping_reason[i] <- "not_stopped"
    }
  }
  
  # Calculate success probability
  summary$success_probability <- mean(success, na.rm = TRUE)
  summary$success_ci <- binom.test(sum(success, na.rm = TRUE), sum(!is.na(success)))$conf.int
  
  # Calculate trial duration statistics
  summary$avg_duration <- mean(trial_duration, na.rm = TRUE)
  summary$median_duration <- median(trial_duration, na.rm = TRUE)
  summary$min_duration <- min(trial_duration, na.rm = TRUE)
  summary$max_duration <- max(trial_duration, na.rm = TRUE)
  summary$sd_duration <- sd(trial_duration, na.rm = TRUE)
  summary$study_durations <- trial_duration
  
  # Calculate hazard ratio statistics
  summary$avg_hr <- mean(hazard_ratios, na.rm = TRUE)
  summary$median_hr <- median(hazard_ratios, na.rm = TRUE)
  summary$min_hr <- min(hazard_ratios, na.rm = TRUE)
  summary$max_hr <- max(hazard_ratios, na.rm = TRUE)
  summary$sd_hr <- sd(hazard_ratios, na.rm = TRUE)
  summary$hazard_ratios <- hazard_ratios
  
  # Calculate power (proportion of p-values less than 0.05)
  summary$power <- mean(p_values < 0.05, na.rm = TRUE)
  
  # Calculate event counts
  summary$avg_events <- mean(event_counts, na.rm = TRUE)
  summary$median_events <- median(event_counts, na.rm = TRUE)
  
  # Calculate early stopping probability
  summary$early_stopping_probability <- mean(early_stopping, na.rm = TRUE)
  
  # Count early stopping reasons
  if (any(early_stopping)) {
    stop_reasons <- table(early_stopping_reason[early_stopping])
    summary$stop_reasons <- as.list(stop_reasons / sum(stop_reasons))
  }
  
  # Create data frame for easy plotting
  summary_df <- data.frame(
    simulation = 1:n_sims,
    success = success,
    duration = trial_duration,
    hazard_ratio = hazard_ratios,
    hr_lower = hr_ci_lower,
    hr_upper = hr_ci_upper,
    p_value = p_values,
    events = event_counts,
    early_stopping = early_stopping,
    stop_reason = early_stopping_reason
  )
  
  summary$summary_df <- summary_df
  
  # Create distribution of treatment effect
  hr_quantiles <- quantile(hazard_ratios, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
  summary$hr_quantiles <- hr_quantiles
  
  # Calculate probability of clinically meaningful effect (HR < 0.8)
  summary$prob_clinical_effect <- mean(hazard_ratios < 0.8, na.rm = TRUE)
  
  # Calculate type I error if truth is no effect (determined by scenario)
  config_scenario <- all_results[[1]]$scenario
  if (config_scenario == "baseline") {
    summary$type_I_error <- summary$success_probability
  }
  
  return(summary)
}

#' @title Plot simulation results
#' 
#' @param sim_result Simulation results object or list of results
#' @param plot_type Type of plot to generate
#' @return ggplot object with the requested visualization
plot_simulation_results <- function(
    sim_result,
    plot_type = c("km", "enrollment", "events", "bayesian")
) {
  plot_type <- match.arg(plot_type)
  
  # Placeholder function - to be implemented
  # This function will generate different types of plots based on simulation results
  
  # This is a placeholder and will be expanded in the future
  if (plot_type == "km") {
    # Return a placeholder plot for now
    p <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, label = "Kaplan-Meier plot will be implemented") +
      theme_minimal() +
      labs(title = "Placeholder for Kaplan-Meier Plot")
  } else {
    # Return a placeholder plot for now
    p <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, label = paste(plot_type, "plot will be implemented")) +
      theme_minimal() +
      labs(title = paste("Placeholder for", plot_type, "Plot"))
  }
  
  return(p)
}

#' @title Export simulation results to various formats
#' 
#' @param sim_result Simulation results object
#' @param format Output format ("csv", "rds", "report")
#' @param file_path Path to save the output
#' @param include_plots Whether to include plots in the report
#' @return Invisible, the path where results were saved
export_simulation_results <- function(
    sim_result,
    format = c("csv", "rds", "report"),
    file_path = NULL,
    include_plots = TRUE
) {
  format <- match.arg(format)
  
  # Generate default file path if not provided
  if (is.null(file_path)) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    
    if (inherits(sim_result, "trial_simulation")) {
      # Single simulation
      trial_name <- sim_result$config$trial_name
      file_name <- paste0(trial_name, "_", timestamp)
    } else {
      # Summary of multiple simulations
      file_name <- paste0("simulation_summary_", timestamp)
    }
    
    base_dir <- "simulation_results"
    if (!dir.exists(base_dir)) {
      dir.create(base_dir, recursive = TRUE)
    }
    
    file_path <- file.path(base_dir, file_name)
  }
  
  if (format == "csv") {
    # Export to CSV format
    if (inherits(sim_result, "trial_simulation")) {
      # Single simulation - export participant data
      write.csv(sim_result$participants, paste0(file_path, "_participants.csv"), row.names = FALSE)
      
      # Export events data
      write.csv(sim_result$events, paste0(file_path, "_events.csv"), row.names = FALSE)
      
      # Export final analysis results
      if (!is.null(sim_result$final_analysis) && is.list(sim_result$final_analysis)) {
        final_analysis_df <- as.data.frame(t(unlist(
          sim_result$final_analysis[!sapply(sim_result$final_analysis, is.list)]
        )))
        write.csv(final_analysis_df, paste0(file_path, "_final_analysis.csv"), row.names = FALSE)
      }
      
      cat("Exported participant data, events, and final analysis to CSV files\n")
    } else {
      # Summary of multiple simulations
      write.csv(sim_result$summary_df, paste0(file_path, ".csv"), row.names = FALSE)
      cat("Exported simulation summary to CSV file\n")
    }
  } else if (format == "rds") {
    # Export to RDS format (R object)
    saveRDS(sim_result, paste0(file_path, ".rds"))
    cat("Exported simulation results to RDS file\n")
  } else if (format == "report") {
    # Create an R Markdown report
    
    # Check if rmarkdown package is available
    if (!requireNamespace("rmarkdown", quietly = TRUE)) {
      stop("The rmarkdown package is required to create reports. Please install it.")
    }
    
    # Create the report content
    report_content <- create_report_content(sim_result, include_plots)
    
    # Write the report to a temporary Rmd file
    rmd_file <- paste0(file_path, ".Rmd")
    writeLines(report_content, rmd_file)
    
    # Render the report to HTML
    html_file <- paste0(file_path, ".html")
    rmarkdown::render(rmd_file, output_file = html_file, quiet = TRUE)
    
    cat("Exported simulation report to HTML file:", html_file, "\n")
  }
  
  return(invisible(file_path))
}

#' @title Create R Markdown report content for simulation results
#' 
#' @param sim_result Simulation results object
#' @param include_plots Whether to include plots
#' @return Character vector with R Markdown content
create_report_content <- function(sim_result, include_plots = TRUE) {
  # Start the R Markdown document
  report <- c(
    "---",
    "title: \"Clinical Trial Simulation Report\"",
    paste0("date: \"", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\""),
    "output: html_document",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
    "library(ggplot2)",
    "library(knitr)",
    "library(survival)",
    "```",
    ""
  )
  
  if (inherits(sim_result, "trial_simulation")) {
    # Single simulation report
    report <- c(report, create_single_sim_report(sim_result, include_plots))
  } else {
    # Multiple simulations summary report
    report <- c(report, create_multi_sim_report(sim_result, include_plots))
  }
  
  return(report)
}

#' @title Create report content for a single simulation
#' 
#' @param sim_result Single simulation result
#' @param include_plots Whether to include plots
#' @return Character vector with R Markdown content
create_single_sim_report <- function(sim_result, include_plots = TRUE) {
  # Extract configuration
  config <- sim_result$config
  
  # Create report sections
  report <- c(
    "## Trial Configuration",
    "",
    paste0("- **Trial Name**: ", config$trial_name),
    paste0("- **Number of Subjects**: ", config$n_subjects),
    paste0("- **Enrollment Rate**: ", config$enrollment_rate, " subjects per time unit"),
    paste0("- **Follow-up Duration**: ", config$follow_up_duration, " time units"),
    paste0("- **Events Required**: ", config$events_required),
    "",
    "### Treatment Arms",
    ""
  )
  
  # Add arm information
  for (arm_name in names(config$arms)) {
    arm <- config$arms[[arm_name]]
    report <- c(report,
      paste0("#### ", arm$name),
      paste0("- **Allocation**: ", arm$allocation * 100, "%"),
      paste0("- **Distribution**: ", arm$tte_distribution),
      "",
      "Parameters:",
      ""
    )
    
    for (param_name in names(arm$tte_params)) {
      report <- c(report,
        paste0("- **", param_name, "**: ", arm$tte_params[[param_name]])
      )
    }
    
    report <- c(report, "")
  }
  
  # Add trial results
  report <- c(report,
    "## Trial Results",
    "",
    paste0("- **Scenario**: ", sim_result$scenario),
    paste0("- **Trial Completion Time**: ", round(sim_result$trial_completion_time, 2), " time units"),
    paste0("- **Event Count**: ", sim_result$event_count),
    "",
    "### Final Analysis",
    ""
  )
  
  # Add final analysis results
  if (!is.null(sim_result$final_analysis) && is.list(sim_result$final_analysis)) {
    # Extract key metrics
    if (!is.null(sim_result$final_analysis$hazard_ratio)) {
      hr <- sim_result$final_analysis$hazard_ratio
      hr_ci_lower <- sim_result$final_analysis$hr_ci_lower
      hr_ci_upper <- sim_result$final_analysis$hr_ci_upper
      p_value <- sim_result$final_analysis$p_value
      
      report <- c(report,
        paste0("- **Hazard Ratio**: ", round(hr, 2), " (95% CI: ", 
               round(hr_ci_lower, 2), " - ", round(hr_ci_upper, 2), ")"),
        paste0("- **p-value**: ", format.pval(p_value, digits = 3)),
        paste0("- **Trial Result**: ", ifelse(p_value < 0.05, "Success", "Failure")),
        ""
      )
    } else {
      report <- c(report,
        "- **Note**: Final analysis details not available",
        ""
      )
    }
  }
  
  # Add interim analyses information
  if (length(sim_result$analyses) > 0) {
    report <- c(report,
      "### Interim Analyses",
      "",
      "| Analysis | Time | Events | Hazard Ratio | p-value | Decision |",
      "|----------|------|--------|--------------|---------|----------|"
    )
    
    for (i in seq_along(sim_result$analyses)) {
      analysis <- sim_result$analyses[[i]]
      
      # Format row data
      time <- round(analysis$time, 2)
      events <- analysis$event_count
      
      if (!is.null(analysis$hazard_ratio)) {
        hr <- round(analysis$hazard_ratio, 2)
        p <- format.pval(analysis$p_value, digits = 3)
      } else {
        hr <- "N/A"
        p <- "N/A"
      }
      
      if (!is.null(analysis$stop_trial)) {
        decision <- ifelse(analysis$stop_trial, 
                          paste0("Stop (", ifelse(analysis$stop_for_efficacy, "efficacy", "futility"), ")"),
                          "Continue")
      } else {
        decision <- "Continue"
      }
      
      report <- c(report,
        paste0("| ", i, " | ", time, " | ", events, " | ", hr, " | ", p, " | ", decision, " |")
      )
    }
    
    report <- c(report, "")
  }
  
  # Add plots if requested
  if (include_plots) {
    report <- c(report,
      "## Visualizations",
      "",
      "### Kaplan-Meier Curve",
      "",
      "```{r km-plot, fig.width=10, fig.height=6}",
      "# Load visualization functions",
      "source('R/visualization.R')",
      "# Create KM plot",
      "plot_km_curve(sim_result)",
      "```",
      "",
      "### Enrollment Over Time",
      "",
      "```{r enrollment-plot, fig.width=10, fig.height=6}",
      "plot_enrollment(sim_result)",
      "```",
      "",
      "### Events Over Time",
      "",
      "```{r events-plot, fig.width=10, fig.height=6}",
      "plot_events(sim_result)",
      "```"
    )
    
    # Add Bayesian plots if available
    if (!is.null(sim_result$bayesian_updates) && 
        is.list(sim_result$bayesian_updates) &&
        !is.null(sim_result$bayesian_updates$hr_posterior)) {
      
      report <- c(report,
        "",
        "### Bayesian Analysis",
        "",
        "```{r bayesian-hr-plot, fig.width=10, fig.height=6}",
        "plot_bayesian_posterior(sim_result, parameter = 'hr')",
        "```",
        "",
        "```{r bayesian-efficacy-plot, fig.width=10, fig.height=6}",
        "plot_bayesian_posterior(sim_result, parameter = 'efficacy')",
        "```"
      )
    }
  }
  
  return(report)
}

#' @title Create report content for multiple simulations
#' 
#' @param sim_summary Summary of multiple simulations
#' @param include_plots Whether to include plots
#' @return Character vector with R Markdown content
create_multi_sim_report <- function(sim_summary, include_plots = TRUE) {
  # Extract configuration from first simulation
  config <- sim_summary$config
  
  # Create report sections
  report <- c(
    "## Simulation Configuration",
    "",
    paste0("- **Trial Configuration**: ", config$trial_name),
    paste0("- **Number of Simulations**: ", sim_summary$n_sims),
    paste0("- **Number of Subjects per Trial**: ", config$n_subjects),
    "",
    "## Simulation Results",
    "",
    "### Key Metrics",
    "",
    paste0("- **Success Probability**: ", round(sim_summary$success_probability * 100, 1), "% (95% CI: ", 
           round(sim_summary$success_ci[1] * 100, 1), "% - ", 
           round(sim_summary$success_ci[2] * 100, 1), "%)"),
    paste0("- **Power**: ", round(sim_summary$power * 100, 1), "%"),
    paste0("- **Probability of Clinically Meaningful Effect (HR < 0.8)**: ", 
           round(sim_summary$prob_clinical_effect * 100, 1), "%"),
    "",
    "### Trial Duration",
    "",
    paste0("- **Average Duration**: ", round(sim_summary$avg_duration, 2), " time units"),
    paste0("- **Median Duration**: ", round(sim_summary$median_duration, 2), " time units"),
    paste0("- **Min Duration**: ", round(sim_summary$min_duration, 2), " time units"),
    paste0("- **Max Duration**: ", round(sim_summary$max_duration, 2), " time units"),
    "",
    "### Treatment Effect",
    "",
    paste0("- **Average Hazard Ratio**: ", round(sim_summary$avg_hr, 2)),
    paste0("- **Median Hazard Ratio**: ", round(sim_summary$median_hr, 2)),
    paste0("- **Hazard Ratio 95% CI**: (", 
           round(sim_summary$hr_quantiles[1], 2), " - ", 
           round(sim_summary$hr_quantiles[5], 2), ")"),
    ""
  )
  
  # Add early stopping information if applicable
  if (!is.null(sim_summary$early_stopping_probability)) {
    report <- c(report,
      "### Early Stopping",
      "",
      paste0("- **Early Stopping Probability**: ", 
             round(sim_summary$early_stopping_probability * 100, 1), "%"),
      ""
    )
    
    if (!is.null(sim_summary$stop_reasons)) {
      report <- c(report, "#### Stopping Reasons", "", "| Reason | Probability |", "|--------|------------|")
      
      for (reason in names(sim_summary$stop_reasons)) {
        prob <- round(sim_summary$stop_reasons[[reason]] * 100, 1)
        report <- c(report, paste0("| ", reason, " | ", prob, "% |"))
      }
      
      report <- c(report, "")
    }
  }
  
  # Add plots if requested
  if (include_plots) {
    report <- c(report,
      "## Visualizations",
      "",
      "### Distribution of Hazard Ratios",
      "",
      "```{r hr-dist-plot, fig.width=10, fig.height=6}",
      "# Load visualization functions",
      "source('R/visualization.R')",
      "# Create hazard ratio distribution plot",
      "plot_simulation_summary(sim_summary, metric = 'hazard_ratio')",
      "```",
      "",
      "### Distribution of Study Duration",
      "",
      "```{r duration-dist-plot, fig.width=10, fig.height=6}",
      "plot_simulation_summary(sim_summary, metric = 'study_duration')",
      "```",
      "",
      "### Success Probability",
      "",
      "```{r success-plot, fig.width=8, fig.height=6}",
      "plot_simulation_summary(sim_summary, metric = 'success_rate')",
      "```"
    )
    
    # Add early stopping plot if applicable
    if (!is.null(sim_summary$early_stopping_probability)) {
      report <- c(report,
        "",
        "### Early Stopping Probability",
        "",
        "```{r early-stopping-plot, fig.width=8, fig.height=6}",
        "plot_simulation_summary(sim_summary, metric = 'early_stopping')",
        "```"
      )
    }
  }
  
  return(report)
}

