# Preamble ----------------------------------------------------------------
# Load packages.
library(tidyverse)
library(rstan)
library(loo)
library(parallel)

# Set the simulation seed.
set.seed(42)

# Load data and ensemble fit.
if (ind_sim == 1) data <- read_rds(here::here("data", str_c("sim_", file_id, ".rds")))
if (ind_emp == 1) data <- read_rds(here::here("data", str_c("emp_", file_id, ".rds")))
ensemble_fit <- read_rds(here::here("output", str_c("ensemble-fit_", file_id, "_", nmember, ".rds")))

# Run the Meta-Learner ----------------------------------------------------
# Create array of likelihoods with effective sample sizes.
LooPSIS_list <- vector(mode = "list", length = nmember)
cores <- detectCores()
for (k in 1:nmember) {
  # Extract log_lik array from each ensemble draws.
  loglik <- ensemble_fit$ensemble_draws[[k]]$log_lik
  ndraw <- dim(loglik)[1]
  LLmat <- matrix(loglik, nr = ndraw)
  
  # Get relative effective sample size for each array.
  r_eff <- relative_eff(x = exp(LLmat), chain_id = double(ndraw) + 1, cores = cores)
  
  # Apply PSIS via loo to array and save.
  LooPSIS_list[[k]] <- loo.matrix(LLmat, r_eff = r_eff, cores = cores, save_psis = FALSE)
}

# Apply the meta-learner and calculate ensemble weights.
ensemble_weights <- loo_model_weights(
  x = LooPSIS_list, 
  method = "stacking", 
  optim_method = "BFGS", 
  optim_control = list(reltol = 1e-10),
  cores = cores
)

# Append weights to model fit.
ensemble_fit$ensemble_weights <- ensemble_weights
write_rds(ensemble_fit, here::here("output", str_c("ensemble-fit_", file_id, "_", nmember, ".rds")))

