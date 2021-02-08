# Preamble ----------------------------------------------------------------
# Load packages.
library(tidyverse)
library(rstan)
library(loo)

# Load Data and Ensemble Fit ----------------------------------------------
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
ensemble_fit <- read_rds(here::here("Output", str_c("ensemble-fit_vb_", file_name, ".rds")))

# Run the Meta-Learner ----------------------------------------------------
# Extract needed draws.
ensemble_draws <- vector(mode = "list", length = length(ensemble_fit))
for (k in 1:length(ensemble_fit)) {
  ensemble_draws[[k]] <- extract(ensemble_fit[[k]], pars = c("Beta", "log_lik"))
}

# Create array of likelihoods with effective sample sizes.
LooPSIS_list <- vector(mode = "list", length = length(ensemble_fit))
cores <- parallel::detectCores()
for (k in 1:length(ensemble_fit)) {
  # Extract log_lik array from each stanfit object.
  LLarray <- ensemble_draws[[k]]$log_lik
  # Get relative effective sample size for each array.
  r_eff <- relative_eff(x = exp(LLarray), cores = cores)
  # Apply PSIS via loo to array and save.
  LooPSIS_list[[k]] = loo.array(LLarray, r_eff = r_eff, cores = cores, save_psis = FALSE)
}

# Calculate weights.
set.seed(22)
ensemble_weights <- loo_model_weights(
  x = LooPSIS_list, 
  method = "stacking", 
  optim_method = "BFGS", 
  optim_control = list(reltol = 1e-10),
  r_eff_list = r_eff_list, # WHERE IS THIS?
  cores = cores
)

# Save weights.
write_rds(ensemble_weights, here::here("Output", str_c("ensemble-weights_", file_name, ".rds")))

