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
ind_ana <- 0        # Indicates attribute non-attendance.
ind_screen <- 1     # Indicates screening.
ind_ana_screen <- 0 # Indicates attribute non-attendance and screening.

hetero <- 1         # Indicates if pathologies differ by individual 

if (ind_none == 1) file_name <- "none"
if (ind_ana == 1) file_name <- "ana"
if (ind_screen == 1) file_name <- "screen"
if (ind_ana_screen == 1) file_name <- "ana-screen"
if (hetero == 1) file_name <- paste(file_name,"-hetero", sep="")
if (hetero == 0) file_name <- paste(file_name,"-homo", sep="")

data <- read_rds(here::here("Data", str_c("sim_", file_name, ".rds")))
hmnl_fit <- read_rds(here::here("Output", str_c("hmnl-fit_", file_name, ".rds")))
ensemble_fit <- read_rds(here::here("Output", str_c("ensemble-fit_vb_", file_name, ".rds")))
ensemble_weights <- read_rds(here::here("Output", str_c("ensemble-weights_", file_name, ".rds")))


# Compute Model Fit -------------------------------------------------------
# Extract needed draws.
hmnl_draws <- extract(hmnl_fit, pars = c("Beta", "Gamma", "Omega", "tau"))
ensemble_draws <- vector(mode = "list", length = length(ensemble_fit))
for (k in 1:length(ensemble_fit)) {
  ensemble_draws[[k]] <- extract(ensemble_fit[[k]], pars = c("Beta", "Gamma", "Omega", "tau", "log_lik"))
}

# Compute HMNL predictive fit.
hmnl_pred_fit <- predictive_fit_hmnl(
  hmnl_fit = hmnl_draws, 
  test_X = data$test_X, 
  test_Y = data$test_Y,
  Z=NULL
)

# Create a model comparision data frame.
model_comparison <- tibble(
  Model = "HMNL",
  LOO = loo(hmnl_fit)$elpd_loo,
  "Hit Rate" = hmnl_pred_fit$hit_rate[2],
  "Hit Prob" = hmnl_pred_fit$hit_prob[2],
)

# Compute fit metrics for ensemble
ensemble_pred_fit <- predictive_fit_ensemble(
  ensemble_weights = ensemble_weights, 
  ensemble_fit = ensemble_draws, 
  test_X = data$test_X, 
  test_Y = data$test_Y,
  Z=NULL
)

# Append results to the model comparison data frame.
model_comparison <- model_comparison %>% 
  bind_rows(
    tibble(
      Model = "Ensemble",
      LOO = ensemble_pred_fit$loo_fit$elpd_loo,
      "Hit Rate" = ensemble_pred_fit$hit_rate[2],
      "Hit Prob" = ensemble_pred_fit$hit_prob[2]
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
write_rds(model_comparison, here::here("Figures", "model_fit.rds"))

