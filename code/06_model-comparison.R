# Preamble ----------------------------------------------------------------
# Load packages.
library(tidyverse)
library(cmdstanr)
library(posterior)
library(loo)

options(mc.cores = 4) # loo option.

# Source functions.
source(here::here("code", "src", "predictive_fit_stacking.R"))
# source(here::here("code", "src", "predictive_fit_hmnl.R"))
# source(here::here("code", "src", "predictive_fit_ensemble.R"))
source(here::here("code", "src", "predictive_fit_hmnl_test.R"))
source(here::here("code", "src", "predictive_fit_ensemble_test.R"))

# Set the simulation seed.
set.seed(42)

# Load data, ensemble, and competing model fit.
if (ind_sim == 1) data <- read_rds(here::here("data", str_c("sim_", file_id, ".rds")))
if (ind_emp == 1) data <- read_rds(here::here("data", str_c("emp_", data_id, ".rds")))
data$test_Z <- matrix(rep(1, nrow(data$test_Y)), ncol = 1)

####################
# Problem importing fit objects with read_rds()? Fixed with fit$save_object()?
# But the ensemble_fit object has a list, not a data frame?
####################

if (ind_sim == 1) hmnl_fit <- read_rds(here::here("output", str_c("hmnl-fit_", file_id, ".rds")))
if (ind_emp == 1) hmnl_fit <- read_rds(here::here("output", str_c("hmnl-fit_", data_id, ".rds")))
ensemble_fit <- read_rds(here::here("output", str_c("ensemble-fit_", file_id, "_", nmember, ".rds")))

####################
# Temp: Estimating the upper bounds.
# HMNL.
hmnl_draws <- hmnl_fit$draws(format = "df", variables = c("Beta", "Gamma", "Omega", "tau"))
hmnl_pred_fit <- predictive_fit_hmnl(
  hmnl_draws = hmnl_draws, 
  test_X = data$test_X, 
  test_Y = data$test_Y,
  test_Z = data$test_Z
)

####################
# Problem with loo working with a cmdstanr object?
####################

# Create a model comparison data frame.
model_comparison <- tibble(
  Model = "HMNL",
  # LOO = loo(hmnl_fit)$elpd_loo, 
  LOO = NA,
  "Hit Rate" = hmnl_pred_fit$hit_rate,
  "Hit Prob" = hmnl_pred_fit$hit_prob
)

# Print results.
model_comparison

# Compute fit metrics for ensemble.
# ensemble_draws <- ensemble_fit$ensemble_draws$draws(format = "df", variables = c("Beta", "Gamma", "Omega", "tau"))
# ensemble_pred_fit <- predictive_fit_hmnl(
#   hmnl_draws = ensemble_fit$ensemble_draws, 
#   test_X = data$test_X, 
#   test_Y = data$test_Y,
#   test_Z = data$test_Z
# )
ensemble_pred_fit <- predictive_fit_ensemble(
  indices = c(ind_none, ind_ana, ind_screen, ind_ana_screen, ind_Z),
  # ensemble_weights = ensemble_fit$ensemble_weights,
  ensemble_weights = 1, # Equal weights.
  ensemble_draws = ensemble_fit$ensemble_draws,
  test_X = data$test_X,
  test_Y = data$test_Y,
  test_Z = data$test_Z,
  mat_ana = ensemble_fit$mat_ana,
  mat_screen = ensemble_fit$mat_screen,
  ensemble_fit = ensemble_fit
)

# Append results to the model comparison data frame.
model_comparison <- model_comparison %>% 
  bind_rows(
    tibble(
      Model = "Ensemble",
      # LOO = ensemble_pred_fit$loo_fit$elpd_loo,
      LOO = NA,
      "Hit Rate" = ensemble_pred_fit$hit_rate,
      "Hit Prob" = ensemble_pred_fit$hit_prob
    )
  )

# Print results.
model_comparison
####################







# Try equal weights.
ensemble_fit$ensemble_weights <- rep(NA, length(ensemble_fit$ensemble_draws))
for (k in 1:length(ensemble_fit$ensemble_weights)) {
  ensemble_fit$ensemble_weights[k] <- 1 / length(ensemble_fit$ensemble_weights)
}

# # Try random weights.
# temp_ensemble_weights <- runif(n = length(ensemble_fit$ensemble_draws), min = 0, max = 1)
# for (k in 1:length(ensemble_fit$ensemble_weights)) {
#   ensemble_fit$ensemble_weights[k] <- temp_ensemble_weights[k] / sum(temp_ensemble_weights)
# }

# # Use the first half of the test data for validation data.
# data$validate_Y <- data$test_Y[1:round(dim(data$test_Y)[1]/2),]
# data$validate_X <- data$test_X[1:round(dim(data$test_X)[1]/2),,,]
# data$validate_Z <- as.matrix(data$test_Z[1:round(dim(data$test_Z)[1]/2),])
# data$test_Y <- data$test_Y[(round(dim(data$test_Y)[1]/2) + 1):dim(data$test_Y)[1],]
# data$test_X <- data$test_X[(round(dim(data$test_X)[1]/2) + 1):dim(data$test_X)[1],,,]
# data$test_Z <- as.matrix(data$test_Z[(round(dim(data$test_Z)[1]/2) + 1):dim(data$test_Z)[1],])
# 
# # Produce predictions for each ensemble member using the validation data.
# ensemble_predictions <- vector(mode = "list", length = nmember)
# for (k in 1:nmember) {
#   ensemble_predictions[[k]] <- predictive_fit_stacking(
#     member_draws = ensemble_fit$ensemble_draws[[k]],
#     validate_X = data$validate_X,
#     validate_Z = data$validate_Z
#   )
# }
# 
# # Produce counts or probabilities to calculate the ensemble weights.
# meta_Y <- as.vector(t(data$validate_Y))
# meta_pred_X <- matrix(NA, nrow = length(meta_Y), ncol = nmember)
# meta_prob_X <- array(NA, dim = c(length(meta_Y), max(meta_Y), nmember))
# for (k in 1:nmember) {
#   for (n in 1:length(meta_Y)) {
#     meta_pred_X[n,k] <- ensemble_predictions[[k]]$predicted_Y[n]
#     meta_prob_X[n,,k] <- ensemble_predictions[[k]]$predicted_probs[,n]
#   }
# }
# temp_ensemble_counts <- rep(NA, nmember)
# temp_ensemble_sum_probs <- rep(NA, nmember)
# for(k in 1:nmember) {
#   temp_ensemble_counts[k] <- sum(meta_Y == meta_pred_X[,k])
#   temp_probs <- NULL
#   for (n in 1:length(meta_Y)) {
#     temp_probs <- c(temp_probs, meta_prob_X[n, meta_Y[n], k])
#   }
#   temp_ensemble_sum_probs[k] <- sum(temp_probs)
# }
# 
# # Normalize the counts or probabilities.
# for (k in 1:nmember) {
#   ensemble_fit$ensemble_weights[k] <- temp_ensemble_counts[k] / sum(temp_ensemble_counts)
#   # ensemble_fit$ensemble_weights[k] <- temp_ensemble_sum_probs[k] / sum(temp_ensemble_sum_probs)
# }

# # Restructure validation data and predictions or probabilities for the meta-learner.
# meta_Y <- as.vector(t(data$validate_Y))
# meta_X <- array(NA, dim = c(length(meta_Y), max(meta_Y), nmember))
# # for (n in 1:length(meta_Y)) {
# #   temp_X <- NULL
# #   for (k in 1:nmember) {
# #     temp_X <- cbind(temp_X, as.vector(t(ensemble_predictions[[k]]$predicted_Y))[n])
# #   }
# #   meta_X[n,,] <- matrix(temp_X, nrow = max(meta_Y), byrow = TRUE)
# # }
# for (k in 1:nmember) {
#   for (n in 1:length(meta_Y)) {
#     meta_X[n,,k] <- ensemble_predictions[[k]]$predicted_probs[,n]
#   }
# }
# 
# # Produce weights for each of the choice tasks in the validation data.
# stan_data <- list(
#   N = dim(meta_X)[1], # Number of observations.
#   A = dim(meta_X)[2], # Number of choice alternatives.
#   L = dim(meta_X)[3], # Number of (estimable) attribute levels.
# 
#   Y = meta_Y,         # Vector of observations.
#   X = meta_X          # Matrix of observation-level covariates.
# )
# 
# options(mc.cores = parallel::detectCores())
# rstan_options(auto_write = FALSE)
# 
# meta_fit <- stan(
#   here::here("Code", "src", "mnl.stan"),
#   data = stan_data,
#   seed = 42
# )
# 
# # Extract weights for each ensemble and normalize.
# temp_ensemble_weights <- extract(meta_fit, pars = c("beta"))
# temp_ensemble_weights <- apply(temp_ensemble_weights$beta, 2, mean)
# temp_ensemble_weights <- (temp_ensemble_weights - min(temp_ensemble_weights)) /
#   (max(temp_ensemble_weights) - min(temp_ensemble_weights))
# for (k in 1:nmember) {
#   ensemble_fit$ensemble_weights[k] <- temp_ensemble_weights[k] / sum(temp_ensemble_weights)
# }

# Check the ensemble weights.
sum(ensemble_fit$ensemble_weights)
hist(ensemble_fit$ensemble_weights)

ensemble_fit$ensemble_weights

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
  mat_screen = ensemble_fit$mat_screen,
  ensemble_fit = ensemble_fit
)

# Append results to the model comparison data frame.
model_comparison <- model_comparison %>% 
  bind_rows(
    tibble(
      Model = "Ensemble",
      # LOO = ensemble_pred_fit$loo_fit$elpd_loo,
      LOO = NA,
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

# # Save model comparison data frame.
# write_rds(model_comparison, here::here("Figures", str_c("model-fit_", file_id, "_", nmember, ".rds")))

# Parameter Recovery ------------------------------------------------------
if (ind_test == 1) {
  # Compare parameter estimates and the true values.
  tibble(
    param = as.factor(1:length(data$bbar)),
    true = data$bbar,
    hmnl_estimate = apply(hmnl_draws$Gamma, 3, mean),
    ensm_estimate_1 = as.vector(ensemble_fit$ensemble_draws[[1]]$Gamma),
    ensm_estimate_2 = as.vector(ensemble_fit$ensemble_draws[[2]]$Gamma)
  ) %>% 
    mutate(
      hmnl_est_diff = true - hmnl_estimate,
      ensm_est_diff_1 = true - ensm_estimate_1,
      ensm_est_diff_2 = true - ensm_estimate_2
    ) %>% 
    pivot_longer(contains("diff"), names_to = "contrast", values_to = "difference") %>%
    ggplot(aes(x = param, y = difference, fill = contrast)) +
    geom_col(position = "dodge") +
    scale_fill_discrete(type = c("gray40", "gray50", "darkred")) +
    # ylim(-2, 2) +
    theme(legend.position = "bottom")
  
  ggsave(here::here("Figures", str_c("parameter-recovery_", file_id, ".png")), width = 7, height = 3)
}

