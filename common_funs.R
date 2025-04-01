# Create survival object
S <- Surv(time, event)  # 0=censored, 1=event

# Fit Kaplan-Meier curve
km_fit <- survfit(Surv(time, event) ~ group, data=df)
plot(km_fit)  # Simple plot
ggsurvplot(km_fit, data=df)  # Enhanced plot with survminer package

# Fit parametric models
weib_model <- survreg(Surv(time, event) ~ covariates, data=df, dist="weibull")
exp_model <- survreg(Surv(time, event) ~ covariates, data=df, dist="exponential")
logn_model <- survreg(Surv(time, event) ~ covariates, data=df, dist="lognormal")

# Extract parameters
weib_shape <- 1/weib_model$scale
weib_scale <- exp(coef(weib_model)[1])

# Generate random survival times from Weibull distribution
# For a clinical trial simulation with covariates
generate_surv_times <- function(n, covariates, model) {
  linear_pred <- model$coefficients[1] + 
    as.matrix(covariates) %*% model$coefficients[-1]
  scale <- exp(linear_pred)
  shape <- 1/model$scale
  rweibull(n, shape=shape, scale=scale)
}

# Simulate censoring times (administrative)
censor_times <- runif(n, min=0, max=study_duration)

# Observed times = minimum of event time and censoring time
observed_times <- pmin(event_times, censor_times)
event_indicators <- as.numeric(event_times <= censor_times)