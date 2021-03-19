# Preamble ----------------------------------------------------------------
# Load packages.
library(tidyverse)
library(rstan)
library(loo)

# Source functions.
source(here::here("Code", "Source", "predictive_fit_hmnl.R"))
source(here::here("Code", "Source", "predictive_fit_ensemble.R"))

# Set the simulation seed.
set.seed(42)

# Load data, ensemble, and competing model fit.
data <- read_rds(here::here("Data", str_c("sim_", file_id, ".rds")))
data$test_Z <- matrix(rep(1, nrow(data$test_Y)), ncol = 1)
hmnl_fit <- read_rds(here::here("Output", str_c("hmnl-fit_", file_id, ".rds")))
ensemble_fit <- read_rds(here::here("Output", str_c("ensemble-fit_", file_id, "_", nmember, ".rds")))
# if (ind_ana == 1) read_rds(here::here("Output", str_c("ana-fit_", file_id, ".rds")))
# if (ind_screen == 1) read_rds(here::here("Output", str_c("screen-fit_", file_id, ".rds")))

# Weed out problems with ELBO for joint-pathology models.
ensemble_fit$ensemble_draws <- 
  ensemble_fit$ensemble_draws[-which(lapply(ensemble_fit$ensemble_draws, length) == 1)]

length(ensemble_fit$ensemble_draws)

# Try equal weights instead.
ensemble_fit$ensemble_weights <- rep(NA, length(ensemble_fit$ensemble_draws)) # For models without any weights.
for (k in 1:length(ensemble_fit$ensemble_weights)) {
  ensemble_fit$ensemble_weights[k] <- 1 / length(ensemble_fit$ensemble_weights)
}

# # Try dropping the final member and renormalizing weights instead.
# ensemble_fit$ensemble_weights[length(ensemble_fit$ensemble_weights)] <- 0
# for (k in 1:length(ensemble_fit$ensemble_weights)) {
#   ensemble_fit$ensemble_weights[k] <- (ensemble_fit$ensemble_weights[k] / sum(ensemble_fit$ensemble_weights))
# }

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

