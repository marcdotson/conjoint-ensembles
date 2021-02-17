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
data <- read_rds(here::here("Data", str_c("sim_", file_id, "_", nmember, ".rds")))
hmnl_fit <- read_rds(here::here("Output", str_c("hmnl-fit_", file_id, ".rds")))
ensemble_fit <- read_rds(here::here("Output", str_c("ensemble-fit_", file_id, "_", nmember, ".rds")))
# if (ind_ana == 1) read_rds(here::here("Output", str_c("ana-fit_", file_id, ".rds")))
# if (ind_screen == 1) read_rds(here::here("Output", str_c("screen-fit_", file_id, ".rds")))

# Compute Model Fit -------------------------------------------------------
# Extract needed draws.
hmnl_draws <- extract(hmnl_fit, pars = c("Beta", "Gamma", "Omega", "tau"))

# Compute HMNL predictive fit.
hmnl_pred_fit <- predictive_fit_hmnl(
  hmnl_fit= hmnl_draws, 
  test_X = data$test_X, 
  test_Y = data$test_Y,
  Z = Z
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
  mat_ana = ensemble_fit$mat_ana,
  mat_screen = ensemble_fit$mat_screen,
  Z=Z
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

#Still Need to add competing models here with predictive fit.

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
# write_rds(model_comparison, here::here("Figures", "model_fit.rds"))

