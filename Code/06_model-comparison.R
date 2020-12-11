# Preamble ----------------------------------------------------------------
# Load packages.
library(tidyverse)
library(rstan)
library(loo)

# Source functions.
source(here::here("Code", "Source", "predictive_fit_ensemble.R"))
# source(here::here("Code", "Source", "hit_rate.R"))
# source(here::here("Code", "Source", "hit_prob.R"))

# Load Data, Fit, and Weights ---------------------------------------------
ind_none <- 0       # Indicates no pathologies.
ind_ana <- 1        # Indicates attribute non-attendance.
ind_screen <- 0     # Indicates screening.
ind_ana_screen <- 0 # Indicates attribute non-attendance and screening.

if (ind_none == 1) file_name <- "none"
if (ind_ana == 1) file_name <- "ana"
if (ind_screen == 1) file_name <- "screen"
if (ind_ana_screen == 1) file_name <- "ana-screen"

data <- read_rds(here::here("Data", str_c("sim_", file_name, ".rds")))
hmnl_fit <- read_rds(here::here("Output", str_c("hmnl-fit_", file_name, ".rds")))
ensemble_fit <- read_rds(here::here("Output", str_c("ensemble-fit_vb_", file_name, ".rds")))
ensemble_weights <- read_rds(here::here("Output", str_c("ensemble-weights_", file_name, ".rds")))

# Compute Model Fit -------------------------------------------------------
# Extract needed draws.
hmnl_draws <- extract(hmnl_fit, pars = "Beta")
ensemble_draws <- vector(mode = "list", length = length(ensemble_fit))
for (k in 1:length(ensemble_fit)) {
  ensemble_draws[[k]] <- extract(ensemble_fit[[k]], pars = c("Beta", "log_lik"))
}

# Compute HMNL predictive fit.
# - Should we have a separate function to generate hold-out sample betas?
# - Do we care about in-sample hit rate and hit prob?

# Create a model comparision data frame.
model_comparison <- tibble(
  Model = "HMNL",
  LOO = loo(hmnl_fit)$elpd_loo,
  "Hit Rate" = NA,
  "Hit Prob" = NA
)

# Compute ensemble predictive fit.
# - Weight log_lik for each ensemble member to get log_lik for ensemble -- in function?
LLarray_ens <- array(0, dim(ensemble_draws[[1]]$log_lik))
cores <- parallel::detectCores()
for (k in 1:length(ensemble_fit)) {
  # Extract log_lik array from each set of draws.
  LLarray_ens <- LLarray_ens + ensemble_weights[k] * ensemble_draws[[k]]$log_lik
}  

# - predictive_fit_ensemble() is using objects from the global environment?
ensemble_pred_fit <- predictive_fit_ensemble(
  ensemble_weights = ensemble_weights, 
  ensemble_fit = ensemble_draws, 
  test_x = data$test_X, 
  test_y = data$test_Y
)

# # Append results to the model comparison data frame.
# model_comparison <- model_comparison %>% 
#   bind_rows(
#     tibble(
#       Model = "Ensemble",
#       LOO = loo_fit_ens$elpd_loo,
#       "Hit Rate" = hit_rate_vec %*% ensemble_weights,
#       "Hit Prob" = hit_prob_vec %*% ensemble_weights
#     )
#   )

# # Save the model comparison data frame.
# write_rds(model_comparison, here::here("Figures", "model_fit.rds"))

