# Preamble ----------------------------------------------------------------
# Load packages.
library(tidyverse)
library(rstan)
library(loo)

# Source functions.
source(here::here("Code", "Source", "predictive_fit_stacking.R"))
source(here::here("Code", "Source", "predictive_fit_hmnl.R"))
source(here::here("Code", "Source", "predictive_fit_ensemble.R"))

# Set the simulation seed.
set.seed(42)

# Load data, ensemble, and competing model fit.
if (ind_sim == 1) data <- read_rds(here::here("Data", str_c("sim_", file_id, ".rds")))
if (ind_emp == 1) data <- read_rds(here::here("Data", str_c("emp_", file_id, ".rds")))
data$test_Z <- matrix(rep(1, nrow(data$test_Y)), ncol = 1)
hmnl_fit <- read_rds(here::here("Output", str_c("hmnl-fit_", file_id, ".rds")))
ensemble_fit <- read_rds(here::here("Output", str_c("ensemble-fit_", file_id, "_", nmember, ".rds")))

# # Try equal weights.
# ensemble_fit$ensemble_weights <- rep(NA, length(ensemble_fit$ensemble_draws))
# for (k in 1:length(ensemble_fit$ensemble_weights)) {
#   ensemble_fit$ensemble_weights[k] <- 1 / length(ensemble_fit$ensemble_weights)
# }

# # Try random weights.
# temp_ensemble_weights <- runif(n = length(ensemble_fit$ensemble_draws), min = 0, max = 1)
# for (k in 1:length(ensemble_fit$ensemble_weights)) {
#   ensemble_fit$ensemble_weights[k] <- temp_ensemble_weights[k] / sum(temp_ensemble_weights)
# }

# Try logistic weights.
# Use the first half of the test data for validation data.
data$validate_Y <- data$test_Y[1:round(dim(data$test_Y)[1]/2),]
data$validate_X <- data$test_X[1:round(dim(data$test_X)[1]/2),,,]
data$validate_Z <- as.matrix(data$test_Z[1:round(dim(data$test_Z)[1]/2),])
data$test_Y <- data$test_Y[(round(dim(data$test_Y)[1]/2) + 1):dim(data$test_Y)[1],]
data$test_X <- data$test_X[(round(dim(data$test_X)[1]/2) + 1):dim(data$test_X)[1],,,]
data$test_Z <- as.matrix(data$test_Z[(round(dim(data$test_Z)[1]/2) + 1):dim(data$test_Z)[1],])

# Produce predictions for each ensemble member using the validation data.
ensemble_predictions <- vector(mode = "list", length = nmember)
for (k in 1:nmember) {
  ensemble_predictions[[k]] <- predictive_fit_stacking(
    member_draws = ensemble_fit$ensemble_draws[[k]],
    validate_X = data$validate_X,
    validate_Z = data$validate_Z
  )
}

# Restructure validation data and predictions for the meta-learner.
meta_Y <- as.vector(t(data$validate_Y))
meta_X <- array(NA, dim = c(length(meta_Y), max(meta_Y), nmember))
for (n in 1:length(meta_Y)) {
  temp_X <- NULL
  for (k in 1:nmember) {
    temp_X <- cbind(temp_X, as.vector(t(ensemble_predictions[[k]]))[n])
  }
  meta_X[n,,] <- matrix(rep(temp_X, max(meta_Y)), nrow = max(meta_Y), byrow = TRUE)
}

# Produce weights for each of the choice tasks in the validation data.
stan_data <- list(
  N = dim(meta_X)[1], # Number of observations.
  A = dim(meta_X)[2], # Number of choice alternatives.
  L = dim(meta_X)[3], # Number of (estimable) attribute levels.
  
  Y = meta_Y,         # Vector of observations.
  X = meta_X          # Matrix of observation-level covariates.
)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = FALSE)

meta_fit <- stan(
  here::here("Code", "Source", "mnl.stan"),
  data = stan_data,
  seed = 42
)

# Extract weights for each ensemble and normalize.
temp_ensemble_weights <- extract(meta_fit, pars = c("beta"))
temp_ensemble_weights <- apply(temp_ensemble_weights$beta, 2, mean)
for (k in 1:nmember) {
  ensemble_fit$ensemble_weights[k] <- temp_ensemble_weights[k] / sum(temp_ensemble_weights)
}

# # Try dropping the final member and normalizing weights.
# ensemble_fit$ensemble_weights[length(ensemble_fit$ensemble_weights)] <- 0
# for (k in 1:length(ensemble_fit$ensemble_weights)) {
#   ensemble_fit$ensemble_weights[k] <- (ensemble_fit$ensemble_weights[k] / sum(ensemble_fit$ensemble_weights))
# }

# # Try weeding out ensemble members with VB errors.
# ensemble_fit$ensemble_draws <-
#   ensemble_fit$ensemble_draws[-which(lapply(ensemble_fit$ensemble_draws, length) == 1)]
# length(ensemble_fit$ensemble_draws)

# Compute Model Fit -------------------------------------------------------
# Extract needed draws.
hmnl_draws <- extract(hmnl_fit, pars = c("Beta", "Gamma", "Omega", "tau"))

# Compute HMNL predictive fit.
hmnl_pred_fit <- predictive_fit_hmnl(
  hmnl_draws = hmnl_draws, 
  test_X = data$test_X, 
  test_Y = data$test_Y,
  test_Z = data$test_Z
)

# Create a model comparison data frame.
model_comparison <- tibble(
  Model = "HMNL",
  LOO = loo(hmnl_fit)$elpd_loo,
  "Hit Rate" = hmnl_pred_fit$hit_rate,
  "Hit Prob" = hmnl_pred_fit$hit_prob
)

# Print results.
model_comparison

# Compute fit metrics for ensemble.
ensemble_pred_fit <- predictive_fit_ensemble(
  indices = c(ind_none, ind_ana, ind_screen, ind_ana_screen, ind_Z),
  ensemble_weights = ensemble_fit$ensemble_weights, 
  ensemble_draws = ensemble_fit$ensemble_draws, 
  test_X = data$test_X, 
  test_Y = data$test_Y,
  test_Z = data$test_Z,
  mat_ana = ensemble_fit$mat_ana,
  mat_screen = ensemble_fit$mat_screen
)

# Append results to the model comparison data frame.
model_comparison <- model_comparison %>% 
  bind_rows(
    tibble(
      Model = "Ensemble",
      LOO = ensemble_pred_fit$loo_fit$elpd_loo,
      "Hit Rate" = ensemble_pred_fit$hit_rate,
      "Hit Prob" = ensemble_pred_fit$hit_prob
    )
  )

# Print results.
model_comparison

# Still need to add competing models here with predictive fit.

# Compute fit metrics for ANA model
#ana_pred_fit <- predictive_fit_ana(
#  ana_fit = ana_draws, 
#  test_X = data$test_X, 
#  test_Y = data$test_Y,
#  Z=NULL
#)

# Append results to the model comparison data frame.
model_comparison <- model_comparison %>% 
  bind_rows(
    tibble(
      Model = "ANA",
      LOO = NA,
      "Hit Rate" = NA,
      "Hit Prob" = NA
    )
  )

#ConjScreen
# Compute fit metrics for Conj Screening model
#cscreen_pred_fit <- predictive_fit_cscreen(
#  cscreen_fit = cscreen_draws, 
#  test_X = data$test_X, 
#  test_Y = data$test_Y,
#  Z=NULL
#)

# Append results to the model comparison data frame.
model_comparison <- model_comparison %>% 
  bind_rows(
    tibble(
      Model = "Conj Screen",
      LOO = NA,
      "Hit Rate" = NA,
      "Hit Prob" = NA
    )
  )

# Save model comparison data frame.
write_rds(model_comparison, here::here("Figures", str_c("model-fit_", file_id, "_", nmember, ".rds")))

