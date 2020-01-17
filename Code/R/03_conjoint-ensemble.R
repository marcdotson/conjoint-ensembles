# Preamble ----------------------------------------------------------------
# Load packages.
library(tidyverse)
library(rstan)
library(loo)
# library(bayesplot)
# library(tidybayes)

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

# Run Ensemble ------------------------------------------------------------
# Specify the data for calibration in a list.
data <- list(
  N = nrow(Y),      # Number of respondents.
  S = ncol(Y),      # Number of choice tasks per respondent.
  P = dim(X)[3],    # Number of product alternatives per choice task.
  L = dim(X)[4],    # Number of (estimable) attribute levels.
  C = ncol(Z),      # Number of respondent-level covariates.
  
  Theta_mean = 0,   # Mean of coefficients for the heterogeneity model.
  Theta_scale = 10, # Scale of coefficients for the heterogeneity model.
  tau_mean = 0,     # Mean of scale parameters for the heterogeneity model.
  tau_scale = 2.5,  # Scale of scale parameters for the heterogeneity model.
  Omega_shape = 2,  # Shape of correlation matrix for the heterogeneity model.
  
  Y = Y,            # Matrix of observed choices.
  X = X,            # Array of experimental designs per choice task.
  Z = Z             # Matrix of respondent-level covariates.
)

# N <- 5; S <- 2
# ns = 0
# for (n in 1:N) {
#   for (s in 1:S) {
#     ns = ns + 1
#     print(ns)
#   }
# }

# Calibrate the model.
fit01 <- stan(
  file = here::here("Code", "Stan", "hmnl_centered.stan"),
  data = data,
  control = list(adapt_delta = 0.99),
  seed = 42
)

# Manual extraction.
log_lik <- extract_log_lik(fit01, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik))
loo(log_lik, r_eff = r_eff)

str(log_lik)
min(log_lik); max(log_lik)

# Fit.
loo(fit01, save_psis = TRUE)

# # Save ensemble output.
# write_rds(fit01, here::here("Output", "hmnl-centered_fit.RDS"))

# Generate Consensus ------------------------------------------------------
