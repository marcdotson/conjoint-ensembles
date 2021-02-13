# Preamble ----------------------------------------------------------------
# Load packages.
library(tidyverse)
library(rstan)
library(loo)

# Source functions.
source(here::here("Code", "Source", "predictive_fit_ensemble.R"))
source(here::here("Code", "Source", "predictive_fit_hmnl.R"))
# source(here::here("Code", "Source", "predictive_fit_ana.R"))
# source(here::here("Code", "Source", "predictive_fit_cscreen.R"))

# Load Data, Fit, and Weights ---------------------------------------------
ind_none <- 0       # Indicates no pathologies.
ind_ana <- 1        # Indicates attribute non-attendance.
ind_screen <- 0     # Indicates screening.
ind_ana_screen <- 0 # Indicates attribute non-attendance and screening.
ind_Z <- 0          # Indicates presence of covariates
nmember <- 400      # Indicate the number of ensemble members.


if (ind_none == 1) file_name <- "none"
if (ind_ana == 1) file_name <- "ana"
if (ind_screen == 1) file_name <- "screen"
if (ind_ana_screen == 1) file_name <- "ana-screen"
if (ind_Z == 1) Z <- data$Z else Z <- NULL

data <- read_rds(here::here("Data", str_c("sim_", file_name, "_", nmember, ".rds")))
hmnl_fit <- read_rds(here::here("Output", str_c("hmnl-fit_", file_name, "_", nmember, ".rds")))
ensemble_draws <- read_rds(here::here("Output", str_c("ensemble-draws_vb_", file_name, "_", nmember, ".rds")))
ensemble_weights <- read_rds(here::here("Output", str_c("ensemble-weights_", file_name, "_", nmember, ".rds")))


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

# Create a model comparision data frame.
model_comparison <- tibble(
  Model = "HMNL",
  LOO = loo(hmnl_fit)$elpd_loo,
  "Hit Rate Gamma Draws" = hmnl_pred_fit$hit_rate_gammadraws,
  "Hit Prob Gamma Draws" = hmnl_pred_fit$hit_prob_gammadraws,
  "Hit Rate Mean of Gammas" = hmnl_pred_fit$hit_rate_meangammas,
  "Hit Prob Mean of Gammas" = hmnl_pred_fit$hit_prob_meangammas
)

# Print results.
model_comparison

# Compute fit metrics for ensemble
ensemble_pred_fit <- predictive_fit_ensemble(
  indices = c(ind_none, ind_ana, ind_screen, ind_ana_screen, ind_Z),
  ensemble_weights = ensemble_weights, 
  ensemble_draws = ensemble_draws, 
  test_X = data$test_X, 
  test_Y = data$test_Y,
  mat_ana = data$mat_ana,
  mat_screen = data$mat_screen,
  Z=Z
)

# Append results to the model comparison data frame.
model_comparison <- model_comparison %>% 
  bind_rows(
    tibble(
      Model = "Ensemble",
      LOO = ensemble_pred_fit$loo_fit$elpd_loo,
      "Hit Rate Gamma Draws" = ensemble_pred_fit$hit_rate_gammadraws,
      "Hit Prob Gamma Draws" = ensemble_pred_fit$hit_prob_gammadraws,
      "Hit Rate Mean of Gammas" = ensemble_pred_fit$hit_rate_meangammas,
      "Hit Prob Mean of Gammas" = ensemble_pred_fit$hit_prob_meangammas
    )
  )

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
      "Hit Rate Gamma Draws" = NA,
      "Hit Prob Gamma Draws" = NA,
      "Hit Rate Mean of Gammas" = NA,
      "Hit Prob Mean of Gammas" = NA
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
      "Hit Rate Gamma Draws" = NA,
      "Hit Prob Gamma Draws" = NA,
      "Hit Rate Mean of Gammas" = NA,
      "Hit Prob Mean of Gammas" = NA
    )
  )

# Save model comparison data frame.
# write_rds(model_comparison, here::here("Figures", str_c("model_fit_", file_name, "_", nmember, ".rds")))
write_rds(model_comparison, here::here("Figures", "model_fit.rds"))
