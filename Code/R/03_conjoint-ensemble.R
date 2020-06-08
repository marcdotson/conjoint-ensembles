# Preamble ----------------------------------------------------------------
# Load packages.
library(tidyverse)
library(rstan)
library(loo)
library(bayesplot)
library(tidybayes)

# Stan options.
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

ind_non <- 1 # Indicates no pathologies.
ind_scr <- 0 # Indicates screening.
ind_ana <- 0 # Indicates attribute non-attendance.
ind_sna <- 0 # Indicates screening and attribute non-attendance.

# Load and format data.
if (ind_non == 1) {
  # Matrix of observed choices.
  Y <- read_csv(here::here("Data", "01_PathologyNone", "Y.csv")) %>% 
    select(-resp) %>% 
    as.matrix()
  
  # Array of experimental designs per choice task.
  X_raw <- read_csv(here::here("Data", "01_PathologyNone", "X.csv"))
}

# Restructure X_raw into X.
X <- array(NA, dim = c(max(X_raw$resp), max(X_raw$task), max(X_raw$alt), ncol(X_raw) - 3))
for (n in 1:max(X_raw$resp)) {
  for (t in 1:max(X_raw$task)) {
    X[n,t,,] <- X_raw %>% 
      filter(resp == n, task == t) %>% 
      select(contains("l_")) %>% 
      as.matrix()
  }
}

# Matrix of respondent-level covariates.
Z <- rep(1, nrow(Y)) %>% 
  as.matrix()

# Randomization -----------------------------------------------------------

# APPLYING CONSTRAINTS DIRECTLY TO PARAMETERS IN STAN?
# NEEDS TO ALLOW FOR DEALING WITH TWO OR MORE PATHOLOGIES JOINTLY.

# Run Ensemble ------------------------------------------------------------
# Calibrate the model with the centered parameterization.
data <- list(
  N = nrow(Y),        # Number of respondents.
  S = ncol(Y),        # Number of choice tasks per respondent.
  P = dim(X)[3],      # Number of product alternatives per choice task.
  L = dim(X)[4],      # Number of (estimable) attribute levels.
  C = ncol(Z),        # Number of respondent-level covariates.
  
  Theta_mean = 0,     # Mean of coefficients for the heterogeneity model.
  Theta_scale = 10,   # Scale of coefficients for the heterogeneity model.
  tau_mean = 0,       # Mean of scale parameters for the heterogeneity model.
  tau_scale = 2.5,    # Scale of scale parameters for the heterogeneity model.
  Omega_shape = 2,    # Shape of correlation matrix for the heterogeneity model.
  
  Y = Y,              # Matrix of observed choices.
  X = X,              # Array of experimental designs per choice task.
  Z = Z               # Matrix of respondent-level covariates.
)

hmnl_centered <- stan(
  file = here::here("Code", "Stan", "hmnl_centered.stan"),
  data = data,
  control = list(adapt_delta = 0.99),
  seed = 42
)

# Save model output.
write_rds(
  hmnl_centered,
  path = here::here("Output", "hmnl_centered.rds")
)

# Calibrate the model with the non-centered parameterization.
data <- list(
  N = nrow(Y),        # Number of respondents.
  S = ncol(Y),        # Number of choice tasks per respondent.
  P = dim(X)[3],      # Number of product alternatives per choice task.
  L = dim(X)[4],      # Number of (estimable) attribute levels.
  C = ncol(Z),        # Number of respondent-level covariates.
  
  Theta_mean = 0,     # Mean of coefficients for the heterogeneity model.
  Theta_scale = 1,    # Scale of coefficients for the heterogeneity model.
  alpha_mean = 0,     # Mean of scale for the heterogeneity model.
  alpha_scale = 10,   # Scale of scale for the heterogeneity model.
  lkj_corr_shape = 5, # Shape of correlation matrix for the heterogeneity model.
  
  Y = Y,              # Matrix of observed choices.
  X = X,              # Array of experimental designs per choice task.
  Z = Z               # Matrix of respondent-level covariates.
)

hmnl_noncentered <- stan(
  file = here::here("Code", "Stan", "hmnl_noncentered.stan"),
  data = data,
  control = list(adapt_delta = 0.99),
  seed = 42
)

# Save model output.
write_rds(
  hmnl_noncentered,
  path = here::here("Output", "hmnl_noncentered.rds")
)

# Generate Consensus ------------------------------------------------------

# CENTERED: 3244.24 seconds (Total)
# NON-CENTERED: 1301.17 seconds (Total)

# Load model output.
hmnl_centered <- read_rds(here::here("Output", "hmnl_centered.rds"))
hmnl_noncentered <- read_rds(here::here("Output", "hmnl_noncentered.rds"))

# Centered parameterization.
log_lik_centered <- extract_log_lik(hmnl_centered, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik_centered))
loo(log_lik_centered, r_eff = r_eff)
loo(hmnl_centered, save_psis = TRUE)

loo_centered <- loo(hmnl_centered, save_psis = TRUE)

# Non-centered parameterization.
log_lik_noncentered <- extract_log_lik(hmnl_noncentered, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik_noncentered))
loo(log_lik_noncentered, r_eff = r_eff)
loo(hmnl_noncentered, save_psis = TRUE)

loo_noncentered <- loo(hmnl_noncentered, save_psis = TRUE)

# Compare model fit.
loo_compare(loo_centered, loo_noncentered)

