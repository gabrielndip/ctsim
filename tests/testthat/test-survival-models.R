context("Survival models")

# Create a test fixture for survival model tests
create_test_survival_data <- function(n = 100) {
  # Create a simple dataset with survival data
  set.seed(123)
  
  data <- data.frame(
    id = 1:n,
    arm = sample(c("control", "treatment"), n, replace = TRUE),
    observed_time = pmax(0.1, rnorm(n, mean = ifelse(rep(1:n, 1) %% 2 == 0, 15, 10), sd = 5)),
    event_status = rbinom(n, 1, 0.7),
    stringsAsFactors = FALSE
  )
  
  return(data)
}

# Test fit_weibull_model function
test_that("fit_weibull_model correctly fits Weibull models to survival data", {
  # Skip if flexsurv package is not available
  skip_if_not_installed("flexsurv")
  
  # Create test data
  data <- create_test_survival_data()
  
  # Fit a Weibull model
  model_fit <- fit_weibull_model(data)
  
  # Check structure
  expect_true(!is.null(model_fit$model))
  expect_true(!is.null(model_fit$params))
  
  # Check that parameters were extracted
  expect_true(!is.null(model_fit$params$shape))
  expect_true(!is.null(model_fit$params$scale_control))
  expect_true(!is.null(model_fit$params$scale_treatment))
  expect_true(!is.null(model_fit$params$hazard_ratio))
  
  # Test with stratification
  data$risk_group <- sample(c("high", "low"), nrow(data), replace = TRUE)
  model_strat <- fit_weibull_model(data, strata = "risk_group")
  
  # Stratified model should have additional terms
  expect_true(length(model_strat$model$coefficients) > length(model_fit$model$coefficients))
})

# Test fit_exponential_model function
test_that("fit_exponential_model correctly fits exponential models to survival data", {
  # Skip if flexsurv package is not available
  skip_if_not_installed("flexsurv")
  
  # Create test data
  data <- create_test_survival_data()
  
  # Fit an exponential model
  model_fit <- fit_exponential_model(data)
  
  # Check structure
  expect_true(!is.null(model_fit$model))
  expect_true(!is.null(model_fit$params))
  
  # Check that parameters were extracted
  expect_true(!is.null(model_fit$params$rate_control))
  expect_true(!is.null(model_fit$params$rate_treatment))
  expect_true(!is.null(model_fit$params$hazard_ratio))
  
  # For exponential, rates should be positive
  expect_true(model_fit$params$rate_control > 0)
  expect_true(model_fit$params$rate_treatment > 0)
})

# Test fit_lognormal_model function
test_that("fit_lognormal_model correctly fits log-normal models to survival data", {
  # Skip if flexsurv package is not available
  skip_if_not_installed("flexsurv")
  
  # Create test data
  data <- create_test_survival_data()
  
  # Fit a log-normal model
  model_fit <- fit_lognormal_model(data)
  
  # Check structure
  expect_true(!is.null(model_fit$model))
  expect_true(!is.null(model_fit$params))
  
  # Check that parameters were extracted
  expect_true(!is.null(model_fit$params$meanlog_control))
  expect_true(!is.null(model_fit$params$meanlog_treatment))
  expect_true(!is.null(model_fit$params$sdlog))
  expect_true(!is.null(model_fit$params$hr_at_median))
})

# Test fit_cox_model function
test_that("fit_cox_model correctly fits Cox proportional hazards models", {
  # Skip if survival package is not available
  skip_if_not_installed("survival")
  
  # Create test data
  data <- create_test_survival_data()
  
  # Fit a Cox model
  model_fit <- fit_cox_model(data)
  
  # Check structure
  expect_true(!is.null(model_fit$model))
  expect_true(!is.null(model_fit$params))
  
  # Check that parameters were extracted
  expect_true(!is.null(model_fit$params$hazard_ratio))
  expect_true(!is.null(model_fit$params$hr_ci_lower))
  expect_true(!is.null(model_fit$params$hr_ci_upper))
  expect_true(!is.null(model_fit$params$p_value))
  
  # Test with stratification
  data$risk_group <- sample(c("high", "low"), nrow(data), replace = TRUE)
  model_strat <- fit_cox_model(data, strata = "risk_group")
  
  # Model should have stratification term
  expect_true(grepl("strata", deparse(model_strat$model$call)))
})

# Test compare_survival_models function
test_that("compare_survival_models compares different survival models", {
  # Skip if required packages are not available
  skip_if_not_installed(c("flexsurv", "survival"))
  
  # Create test data
  data <- create_test_survival_data()
  
  # Compare models
  comparison <- compare_survival_models(data)
  
  # Check structure
  expect_true(!is.null(comparison$models))
  expect_true(!is.null(comparison$metrics))
  expect_true(!is.null(comparison$predictions))
  
  # Should have all model types
  expect_named(comparison$models, c("weibull", "exponential", "lognormal", "cox"))
  
  # Metrics should be a data frame with AIC and BIC
  expect_true(is.data.frame(comparison$metrics))
  expect_true("AIC" %in% names(comparison$metrics))
  expect_true("BIC" %in% names(comparison$metrics))
  
  # Predictions should have control and treatment columns
  expect_true(is.data.frame(comparison$predictions))
  expect_true("control" %in% names(comparison$predictions))
  expect_true("treatment" %in% names(comparison$predictions))
})

# Test generate_survival_times function
test_that("generate_survival_times creates random variates from survival distributions", {
  # Test Weibull distribution
  times_weibull <- generate_survival_times(
    n = 100,
    distribution = "weibull",
    params = list(shape = 1.2, scale = 12)
  )
  
  # Check length and positivity
  expect_length(times_weibull, 100)
  expect_true(all(times_weibull > 0))
  
  # Test exponential distribution
  times_exp <- generate_survival_times(
    n = 100,
    distribution = "exponential",
    params = list(rate = 0.1)
  )
  
  expect_length(times_exp, 100)
  expect_true(all(times_exp > 0))
  
  # Test log-normal distribution
  times_lnorm <- generate_survival_times(
    n = 100,
    distribution = "lognormal",
    params = list(meanlog = 2, sdlog = 0.5)
  )
  
  expect_length(times_lnorm, 100)
  expect_true(all(times_lnorm > 0))
  
  # Test gamma distribution
  times_gamma <- generate_survival_times(
    n = 100,
    distribution = "gamma",
    params = list(shape = 2, rate = 0.1)
  )
  
  expect_length(times_gamma, 100)
  expect_true(all(times_gamma > 0))
  
  # Test error for missing parameters
  expect_error(
    generate_survival_times(n = 10, distribution = "weibull", params = list(shape = 1.2)),
    "Weibull distribution requires shape and scale parameters"
  )
  
  # Test error for invalid distribution
  expect_error(
    generate_survival_times(n = 10, distribution = "invalid"),
    "should be one of"
  )
})

# Test create_competing_risks_model function
test_that("create_competing_risks_model builds valid competing risks generator", {
  # Define competing risks
  cr_model <- create_competing_risks_model(
    cause_distributions = list(
      "primary" = "weibull",
      "secondary" = "exponential"
    ),
    cause_params = list(
      "primary" = list(shape = 1.2, scale = 12),
      "secondary" = list(rate = 1/24)
    )
  )
  
  # Check that function is returned
  expect_type(cr_model, "closure")
  
  # Generate data
  cr_data <- cr_model(n = 50)
  
  # Check structure
  expect_true(is.data.frame(cr_data))
  expect_equal(nrow(cr_data), 50)
  expect_true("time" %in% names(cr_data))
  expect_true("cause" %in% names(cr_data))
  
  # Times should be positive
  expect_true(all(cr_data$time > 0))
  
  # Causes should be either "primary" or "secondary"
  expect_true(all(cr_data$cause %in% c("primary", "secondary")))
  
  # Test error case: mismatched inputs
  expect_error(
    create_competing_risks_model(
      cause_distributions = list("primary" = "weibull"),
      cause_params = list("primary" = list(shape = 1.2, scale = 12), "extra" = list())
    ),
    "cause_distributions and cause_params must have the same length"
  )
})

# Test create_multistate_model function
test_that("create_multistate_model builds valid multi-state generator", {
  # Define a simple multi-state model
  msm_model <- create_multistate_model(
    states = c("Healthy", "Disease", "Death"),
    transitions = matrix(c(
      0, 0.1, 0.05,
      0, 0, 0.2,
      0, 0, 0
    ), nrow = 3, byrow = TRUE),
    arm_effect = list(
      "control" = 1,
      "treatment" = 0.7
    )
  )
  
  # Check that function is returned
  expect_type(msm_model, "closure")
  
  # Generate data for control arm
  msm_data_control <- msm_model(n = 30, arm = "control", max_time = 24)
  
  # Check structure
  expect_true(is.data.frame(msm_data_control))
  expect_true("id" %in% names(msm_data_control))
  expect_true("time" %in% names(msm_data_control))
  expect_true("state" %in% names(msm_data_control))
  
  # Generate data for treatment arm
  msm_data_treatment <- msm_model(n = 30, arm = "treatment", max_time = 24)
  
  # Treatment should generally have fewer transitions due to lower intensity
  # But this is probabilistic, so we can't test it directly without many samples
  
  # Test arm effect with custom arm
  msm_data_custom <- msm_model(n = 30, arm = "custom", max_time = 24)
  # Should use default effect of 1
  expect_true(nrow(msm_data_custom) > 0)
})