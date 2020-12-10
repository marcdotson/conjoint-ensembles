# Load packages.
library(tidyverse)
library(rstan)
# library(bayesplot)
# library(tidybayes)
library(loo)

# Source fit functions.
source(here::here("Code", "Source", "predictive_fit_ensemble.R"))
source(here::here("Code", "Source", "hit_rate.R"))
source(here::here("Code", "Source", "hit_prob.R"))

ind_none <- 0       # Indicates no pathologies.
ind_ana <- 1        # Indicates attribute non-attendance.
ind_screen <- 0     # Indicates screening.
ind_ana_screen <- 0 # Indicates attribute non-attendance and screening.
ind_real <- 0       # Indicates ____ data.

if (ind_none == 1) file_name <- "none"
if (ind_ana == 1) file_name <- "ana"
if (ind_screen == 1) file_name <- "screen"
if (ind_ana_screen == 1) file_name <- "ana-screen"
if (ind_real == 1) file_name <- "design"

# Load hmnl-fit output.
hmnl_fit <- read_rds(here::here("Output", str_c("hmnl-fit_", file_name, ".rds")))

# Compute fit.
hmnl_draws <- extract(hmnl_fit, pars = "Beta")

model_comparison <- tibble(
  Model = "HMNL",
  LOO = loo(hmnl_fit)$elpd_loo,
  "Hit Rate" = hit_rate_vec,
  "Hit Prob" = hit_prob_vec
)

# Save as a data frame with ensemble-fit and competing model fit.

# Load weighted enemble-fit output.
out <- read_rds(here::here("Output", str_c("weighted_ensemble-fit_vb_", file_name, ".rds")))

#prepare test data
test_Y <- out$test_Y
test_X <- out$test_X

# weights <- out$weights
# dims <- dim(ensemble_fit$ensemble_fit[[1]]$log_lik) #dimensions of log_lik values
cores <- 4 # number of cores

#weight log_lik for each model to get log_lik for ensemble
LLarray_ens = array(0,dims)
for(k in 1:n_ens){
  #extract log_lik array from each stanfit object
  LLarray_ens <- LLarray_ens + weights[k]*ensemble_fit[[k]]$log_lik
}  

# #get functions for predictive fit
# ana_fit <- predictive_fit_ensemble(ensemble_weights=weights, ensemble_fit=ensemble_fit, test_X, test_Y)


model_comparison <- model_comparison %>% 
  bind_rows(
    tibble(
      Model = "Ensemble",
      LOO = loo_fit_ens$elpd_loo,
      "Hit Rate" = hit_rate_vec%*%ensemble_weights,
      "Hit Prob" = hit_prob_vec%*%ensemble_weights
    )
  )

write_rds(model_comparison, here::here("Figures", "model_fit.rds"))

