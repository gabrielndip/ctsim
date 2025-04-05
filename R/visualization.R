# Clinical Trial Simulation
# visualization.R
# Created: 2025-03-31
# Description: Visualization functions for clinical trial simulations

library(tidyverse)
library(ggplot2)
library(survival)
library(gridExtra)

#' @title Create Kaplan-Meier survival curve plot
#' 
#' @param sim_result Simulation results object
#' @param conf_int Logical, whether to show confidence intervals
#' @param risk_table Logical, whether to show risk table
#' @param title Plot title
#' @return ggplot object with Kaplan-Meier plot
plot_km_curve <- function(
    sim_result,
    conf_int = TRUE,
    risk_table = TRUE,
    title = "Kaplan-Meier Survival Curves"
) {
  # Extract participant data
  participants <- sim_result$participants
  
  # Fit Kaplan-Meier model
  km_fit <- survfit(Surv(observed_time, event_status) ~ arm, data = participants)
  
  # Create plot data frame
  km_df <- data.frame(
    time = km_fit$time,
    surv = km_fit$surv,
    upper = km_fit$upper,
    lower = km_fit$lower,
    strata = km_fit$strata
  )
  
  # Create color-blind friendly colors
  arm_colors <- c("control" = "#1F78B4", "treatment" = "#33A02C")
  
  # Generate base plot
  p <- ggplot(km_df, aes(x = time, y = surv, color = strata)) +
    geom_step(linewidth = 1) +
    scale_color_manual(
      name = "Study Arm",
      values = arm_colors,
      labels = c("Control", "Treatment")
    ) +
    labs(
      title = title,
      x = "Time",
      y = "Survival Probability"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    ) +
    coord_cartesian(ylim = c(0, 1))
  
  # Add confidence intervals if requested
  if (conf_int) {
    p <- p + geom_ribbon(
      aes(ymin = lower, ymax = upper, fill = strata),
      alpha = 0.2, linetype = 0
    ) +
    scale_fill_manual(
      name = "Study Arm",
      values = arm_colors,
      labels = c("Control", "Treatment")
    )
  }
  
  # Add censoring marks
  cens_data <- data.frame(
    time = km_fit$time[km_fit$n.censor > 0],
    surv = km_fit$surv[km_fit$n.censor > 0],
    n.censor = km_fit$n.censor[km_fit$n.censor > 0],
    strata = km_fit$strata[km_fit$n.censor > 0]
  )
  
  if (nrow(cens_data) > 0) {
    p <- p + geom_point(
      data = cens_data,
      aes(x = time, y = surv),
      shape = 3, size = 2
    )
  }
  
  # Add risk table if requested
  if (risk_table) {
    # Create risk table data
    at_risk <- summary(km_fit, times = seq(0, max(participants$observed_time), length.out = 6))
    risk_data <- data.frame(
      time = rep(at_risk$time, 2),
      arm = rep(c("Control", "Treatment"), each = length(at_risk$time)),
      n.risk = c(at_risk$n.risk[1:length(at_risk$time)], 
                at_risk$n.risk[(length(at_risk$time)+1):(length(at_risk$time)*2)])
    )
    
    # Create risk table plot
    risk_plot <- ggplot(risk_data, aes(x = time, y = arm, label = n.risk)) +
      geom_text(size = 3) +
      theme_minimal() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10, hjust = 1),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
      ) +
      labs(title = "Number at risk") +
      scale_x_continuous(breaks = at_risk$time, labels = round(at_risk$time, 1))
    
    # Combine plots
    grid_plots <- grid.arrange(p, risk_plot, heights = c(4, 1), ncol = 1)
    return(grid_plots)
  }
  
  return(p)
}

#' @title Plot enrollment over time
#' 
#' @param sim_result Simulation results object
#' @param by_arm Logical, whether to stratify by arm
#' @param title Plot title
#' @return ggplot object with enrollment plot
plot_enrollment <- function(
    sim_result,
    by_arm = TRUE,
    title = "Enrollment Over Time"
) {
  # Extract event data
  events <- sim_result$events
  
  # Filter for enrollment events
  enrollment_events <- events[events$event_type == "enrollment", ]
  
  # Sort by time
  enrollment_events <- enrollment_events[order(enrollment_events$time), ]
  
  # Create cumulative enrollment
  enrollment_events$cumulative <- 1:nrow(enrollment_events)
  
  # Create plot
  p <- ggplot(enrollment_events, aes(x = time, y = cumulative)) +
    labs(
      title = title,
      x = "Time",
      y = "Cumulative Enrollment"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  # Add enrollment by arm if requested
  if (by_arm) {
    # Calculate cumulative enrollment by arm
    arms <- unique(enrollment_events$arm)
    for (arm in arms) {
      arm_events <- enrollment_events[enrollment_events$arm == arm, ]
      arm_events <- arm_events[order(arm_events$time), ]
      arm_events$arm_cumulative <- 1:nrow(arm_events)
      enrollment_events$arm_cumulative[enrollment_events$arm == arm] <- arm_events$arm_cumulative
    }
    
    # Use arm colors
    arm_colors <- c("control" = "#1F78B4", "treatment" = "#33A02C")
    
    p <- p + 
      geom_step(color = "gray50", linewidth = 1) +
      geom_step(aes(y = arm_cumulative, color = arm), linewidth = 1) +
      scale_color_manual(
        name = "Study Arm",
        values = arm_colors,
        labels = c("Control", "Treatment")
      )
  } else {
    p <- p + geom_step(color = "#1F78B4", linewidth = 1)
  }
  
  return(p)
}

#' @title Plot events over time
#' 
#' @param sim_result Simulation results object
#' @param by_arm Logical, whether to stratify by arm
#' @param include_censoring Logical, whether to include censoring events
#' @param title Plot title
#' @return ggplot object with events plot
plot_events <- function(
    sim_result,
    by_arm = TRUE,
    include_censoring = FALSE,
    title = "Events Over Time"
) {
  # Extract event data
  events <- sim_result$events
  
  # Filter for relevant events
  if (include_censoring) {
    relevant_events <- events[events$event_type %in% c("event", "censoring"), ]
  } else {
    relevant_events <- events[events$event_type == "event", ]
  }
  
  # Sort by time
  relevant_events <- relevant_events[order(relevant_events$time), ]
  
  # Create cumulative events
  if (include_censoring) {
    # Create separate counts for events and censoring
    event_counts <- relevant_events$event_count
    censor_counts <- relevant_events$censored_count
    relevant_events$cumulative <- event_counts + censor_counts
  } else {
    relevant_events$cumulative <- relevant_events$event_count
  }
  
  # Create plot
  p <- ggplot(relevant_events, aes(x = time, y = cumulative)) +
    labs(
      title = title,
      x = "Time",
      y = "Cumulative Events"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  # Add required events threshold
  if (!is.null(sim_result$config$events_required)) {
    p <- p + geom_hline(
      yintercept = sim_result$config$events_required,
      linetype = "dashed",
      color = "red"
    ) +
    annotate(
      "text",
      x = min(relevant_events$time) + 0.1 * (max(relevant_events$time) - min(relevant_events$time)),
      y = sim_result$config$events_required * 1.05,
      label = paste("Required events:", sim_result$config$events_required),
      color = "red",
      hjust = 0
    )
  }
  
  # Add events by arm if requested
  if (by_arm && !include_censoring) {
    # For simplicity, we'll only show events by arm if censoring is not included
    arm_events <- data.frame()
    
    # Calculate cumulative events by arm
    for (arm_name in names(sim_result$config$arms)) {
      arm_data <- relevant_events[relevant_events$arm == arm_name, ]
      if (nrow(arm_data) > 0) {
        arm_data$arm_name <- arm_name
        arm_data$arm_cumulative <- cumsum(arm_data$event_type == "event")
        arm_events <- rbind(arm_events, arm_data)
      }
    }
    
    if (nrow(arm_events) > 0) {
      # Use arm colors
      arm_colors <- c("control" = "#1F78B4", "treatment" = "#33A02C")
      
      p <- p + 
        geom_step(color = "gray50", linewidth = 1) +
        geom_step(data = arm_events, aes(y = arm_cumulative, color = arm_name), linewidth = 1) +
        scale_color_manual(
          name = "Study Arm",
          values = arm_colors,
          labels = c("Control", "Treatment")
        )
    }
  } else {
    event_color <- ifelse(include_censoring, "gray50", "#1F78B4")
    p <- p + geom_step(color = event_color, linewidth = 1)
    
    # Add separate lines for events and censoring if requested
    if (include_censoring) {
      p <- p + 
        geom_step(aes(y = event_count), color = "#1F78B4", linewidth = 1) +
        geom_step(aes(y = censored_count), color = "#33A02C", linewidth = 1, linetype = "dashed") +
        annotate(
          "text",
          x = max(relevant_events$time) * 0.8,
          y = max(relevant_events$event_count) * 0.5,
          label = "Events",
          color = "#1F78B4"
        ) +
        annotate(
          "text",
          x = max(relevant_events$time) * 0.8,
          y = max(relevant_events$censored_count) * 0.5,
          label = "Censored",
          color = "#33A02C"
        )
    }
  }
  
  return(p)
}

#' @title Plot Bayesian posterior distributions
#' 
#' @param sim_result Simulation results object
#' @param parameter Parameter to visualize ("hr", "efficacy")
#' @param title Plot title
#' @return ggplot object with posterior distribution plot
plot_bayesian_posterior <- function(
    sim_result,
    parameter = c("hr", "efficacy"),
    title = NULL
) {
  parameter <- match.arg(parameter)
  
  # Check if Bayesian updates exist
  if (is.null(sim_result$bayesian_updates) || 
      !is.list(sim_result$bayesian_updates) ||
      is.null(sim_result$bayesian_updates$hr_posterior)) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "No Bayesian analysis available") +
             theme_void())
  }
  
  posterior <- sim_result$bayesian_updates$hr_posterior
  
  if (parameter == "hr") {
    # Plot hazard ratio posterior
    hr_samples <- exp(posterior$samples$b_armtreatment)
    hr_df <- data.frame(hazard_ratio = hr_samples)
    
    # Set default title if not provided
    if (is.null(title)) {
      title <- "Posterior Distribution of Hazard Ratio"
    }
    
    # Create plot
    p <- ggplot(hr_df, aes(x = hazard_ratio)) +
      geom_histogram(fill = "#1F78B4", color = "white", bins = 30) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
      geom_vline(xintercept = exp(posterior$median), linetype = "solid", color = "black") +
      annotate(
        "text",
        x = exp(posterior$median),
        y = 0,
        label = paste("Median HR:", round(exp(posterior$median), 2)),
        vjust = -1
      ) +
      labs(
        title = title,
        x = "Hazard Ratio (Treatment vs Control)",
        y = "Frequency"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5)
      )
  } else if (parameter == "efficacy") {
    # Plot probability of efficacy
    if (is.null(sim_result$bayesian_updates$prob_efficacy)) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, label = "No efficacy probability available") +
               theme_void())
    }
    
    efficacy_prob <- sim_result$bayesian_updates$prob_efficacy
    superiority_prob <- sim_result$bayesian_updates$prob_superiority
    
    # Set default title if not provided
    if (is.null(title)) {
      title <- "Probability of Treatment Efficacy"
    }
    
    # Create data for bar plot
    prob_df <- data.frame(
      probability_type = c("Clinically Meaningful Effect (HR < 0.8)", 
                           "Any Benefit (HR < 1)"),
      probability = c(efficacy_prob, superiority_prob)
    )
    
    # Create plot
    p <- ggplot(prob_df, aes(x = probability_type, y = probability, fill = probability_type)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = paste0(round(probability * 100, 1), "%")), 
                vjust = -0.5, size = 5) +
      scale_fill_manual(values = c("#1F78B4", "#33A02C")) +
      labs(
        title = title,
        x = "",
        y = "Probability"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      ) +
      ylim(0, 1)
  }
  
  return(p)
}

#' @title Plot simulation summary across multiple simulations
#' 
#' @param sim_summary Summary of multiple simulation results
#' @param metric Metric to visualize
#' @param title Plot title
#' @return ggplot object with summary plot
plot_simulation_summary <- function(
    sim_summary,
    metric = c("success_rate", "early_stopping", "study_duration", "hazard_ratio"),
    title = NULL
) {
  metric <- match.arg(metric)
  
  # Default titles
  default_titles <- list(
    success_rate = "Probability of Trial Success",
    early_stopping = "Early Stopping Probability",
    study_duration = "Distribution of Study Duration",
    hazard_ratio = "Distribution of Estimated Hazard Ratios"
  )
  
  # Use default title if not provided
  if (is.null(title)) {
    title <- default_titles[[metric]]
  }
  
  # Create appropriate plot based on metric
  if (metric == "success_rate" || metric == "early_stopping") {
    # For binary outcomes, create probability bar plot
    if (metric == "success_rate") {
      value <- sim_summary$success_probability
      label <- paste0(round(value * 100, 1), "%")
    } else {
      value <- sim_summary$early_stopping_probability
      label <- paste0(round(value * 100, 1), "%")
    }
    
    df <- data.frame(value = value)
    
    p <- ggplot(df, aes(x = 1, y = value)) +
      geom_bar(stat = "identity", fill = "#1F78B4", width = 0.5) +
      geom_text(aes(label = label), vjust = -0.5, size = 5) +
      labs(
        title = title,
        x = "",
        y = "Probability"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      ) +
      ylim(0, 1)
    
  } else if (metric == "study_duration") {
    # For continuous outcomes, create histogram or density plot
    if (is.null(sim_summary$study_durations)) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, label = "No study duration data available") +
               theme_void())
    }
    
    df <- data.frame(duration = sim_summary$study_durations)
    
    p <- ggplot(df, aes(x = duration)) +
      geom_histogram(fill = "#1F78B4", color = "white", bins = 20) +
      geom_vline(xintercept = mean(df$duration), linetype = "dashed", color = "red") +
      annotate(
        "text",
        x = mean(df$duration),
        y = 0,
        label = paste("Mean:", round(mean(df$duration), 1)),
        vjust = -1
      ) +
      labs(
        title = title,
        x = "Study Duration",
        y = "Frequency"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5)
      )
    
  } else if (metric == "hazard_ratio") {
    # For hazard ratio distribution
    if (is.null(sim_summary$hazard_ratios)) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, label = "No hazard ratio data available") +
               theme_void())
    }
    
    df <- data.frame(hazard_ratio = sim_summary$hazard_ratios)
    
    p <- ggplot(df, aes(x = hazard_ratio)) +
      geom_histogram(fill = "#1F78B4", color = "white", bins = 20) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
      geom_vline(xintercept = mean(df$hazard_ratio), linetype = "solid", color = "black") +
      annotate(
        "text",
        x = mean(df$hazard_ratio),
        y = 0,
        label = paste("Mean:", round(mean(df$hazard_ratio), 2)),
        vjust = -1
      ) +
      labs(
        title = title,
        x = "Hazard Ratio (Treatment vs Control)",
        y = "Frequency"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5)
      )
  }
  
  return(p)
}

#' @title Create a dashboard of multiple plots for a simulation
#' 
#' @param sim_result Simulation results object
#' @param plots Vector of plot types to include
#' @param ncol Number of columns in the grid
#' @param title Main dashboard title
#' @return A grid of plots
create_simulation_dashboard <- function(
    sim_result,
    plots = c("km", "enrollment", "events", "bayesian_hr"),
    ncol = 2,
    title = "Clinical Trial Simulation Dashboard"
) {
  # Create list to hold all plots
  plot_list <- list()
  
  # Create each requested plot
  for (plot_type in plots) {
    if (plot_type == "km") {
      p <- plot_km_curve(sim_result, conf_int = TRUE, risk_table = FALSE)
      plot_list <- c(plot_list, list(p))
    } else if (plot_type == "enrollment") {
      p <- plot_enrollment(sim_result, by_arm = TRUE)
      plot_list <- c(plot_list, list(p))
    } else if (plot_type == "events") {
      p <- plot_events(sim_result, by_arm = TRUE)
      plot_list <- c(plot_list, list(p))
    } else if (plot_type == "bayesian_hr") {
      p <- plot_bayesian_posterior(sim_result, parameter = "hr")
      plot_list <- c(plot_list, list(p))
    } else if (plot_type == "bayesian_efficacy") {
      p <- plot_bayesian_posterior(sim_result, parameter = "efficacy")
      plot_list <- c(plot_list, list(p))
    }
  }
  
  # Combine plots into a grid
  nrow <- ceiling(length(plot_list) / ncol)
  
  # Add title
  title_grob <- grid::textGrob(
    title,
    gp = grid::gpar(fontface = "bold", fontsize = 16)
  )
  
  # Arrange plots in grid
  grid_arranged <- gridExtra::grid.arrange(
    grobs = plot_list,
    ncol = ncol,
    nrow = nrow,
    top = title_grob
  )
  
  return(grid_arranged)
}

