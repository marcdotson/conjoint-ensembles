# Preamble ----------------------------------------------------------------
# Load packages.
library(tidyverse)
library(rstan)
library(bayesplot)
library(tidybayes)
library(loo)

# Set Stan options.
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# # Load data.
# ind_non <- 1 # Indicates no pathologies.
# ind_scr <- 0 # Indicates screening.
# ind_ana <- 0 # Indicates attribute non-attendance.
# ind_sna <- 0 # Indicates screening and attribute non-attendance.
# 
# # Load and format data.
# if (ind_non == 1) {
#   # Matrix of observed choices.
#   Y <- read_csv(here::here("Data", "01_PathologyNone", "Y.csv")) %>%
#     select(-resp) %>%
#     as.matrix()
# 
#   # Array of experimental designs per choice task.
#   X_raw <- read_csv(here::here("Data", "01_PathologyNone", "X.csv"))
# }
# 
# # Restructure X_raw into X.
# X <- array(NA, dim = c(max(X_raw$resp), max(X_raw$task), max(X_raw$alt), ncol(X_raw) - 3))
# for (n in 1:max(X_raw$resp)) {
#   for (t in 1:max(X_raw$task)) {
#     X[n,t,,] <- X_raw %>%
#       filter(resp == n, task == t) %>%
#       select(contains("l_")) %>%
#       as.matrix()
#   }
# }
# 
# # Matrix of respondent-level covariates.
# Z <- rep(1, nrow(Y)) %>%
#   as.matrix()

# Jeff's data (not working) but still using the indicator matrices.
# load(here::here("Data", "R05_Zerorez", "design.RData"))
load(here::here("Data", "R06_Ground-Beef", "design.RData"))
Y <- data$train_Y
X <- data$train_X
Z <- matrix(rep(1, nrow(Y)), ncol = 1)

# Work with just two members of the ensemble to start.
ind_ana <- data$mat_ana[1:2, 1:dim(X)[4]]
data$mat_ana <- ind_ana
ind_screen <- data$mat_screen[1:2, 1:dim(X)[4]]
data$mat_screen <- ind_screen

# Initialize Ensemble with Posterior Means --------------------------------
stan_data <- list(
  R = dim(X)[1],    # Number of respondents.
  S = dim(X)[2],    # Number of choice tasks.
  A = dim(X)[3],    # Number of choice alternatives.
  I = dim(X)[4],    # Number of observation-level covariates.
  J = ncol(Z),      # Number of population-level covariates.
  
  Gamma_mean = 0,   # Mean of population-level means.
  Gamma_scale = 5,  # Scale of population-level means.
  Omega_shape = 2,  # Shape of population-level scale.
  tau_mean = 0,     # Mean of population-level scale.
  tau_scale = 5,    # Scale of population-level scale.
  
  Y = Y,            # Matrix of observations.
  X = X,            # Array of observation-level covariates.
  Z = Z             # Matrix of population-level covariates.
)

initial_fit <- stan(
  here::here("Code", "Source", "hmnl.stan"),
  data = stan_data,
  seed = 42
)

# Save initial fit output.
write_rds(
  initial_fit,
  here::here("Output", "initial_fit_pathology-none.rds")
)

initial_values <- extract(initial_fit, pars = c("Gamma", "Omega", "tau", "Delta", "Beta"))

# Set initial values by providing a list equal in length to the number of chains.
# The elements of this list should themselves be named lists, where each of
# these named lists has the name of a parameter and is used to specify the
# initial values for that parameter for the corresponding chain.

# Ensemble Calibration ----------------------------------------------------
K <- nrow(ind_ana)
ensemble_fit <- vector(mode = "list", length = K)
for (k in 1:K) {
  stan_data <- list(
    R = dim(X)[1],    # Number of respondents.
    S = dim(X)[2],    # Number of choice tasks.
    A = dim(X)[3],    # Number of choice alternatives.
    I = dim(X)[4],    # Number of observation-level covariates.
    J = ncol(Z),      # Number of population-level covariates.
    K = K,            # Number of members of the ensemble.
    k = k,            # Ensemble member number.
    
    Gamma_mean = 0,   # Mean of population-level means.
    Gamma_scale = 5,  # Scale of population-level means.
    Omega_shape = 2,  # Shape of population-level scale.
    tau_mean = 0,     # Mean of population-level scale.
    tau_scale = 5,    # Scale of population-level scale.
    
    Y = Y,            # Matrix of observations.
    X = X,            # Array of observation-level covariates.
    Z = Z,            # Matrix of population-level covariates.
    ind_ana = ind_ana # Matrix of ensemble indicators for ANA.
  )
  
  ensemble_fit[[k]] <- vb(
    stan_model(here::here("Code", "Source", "hmnl_ensemble.stan")),
    data = stan_data,
    init = initial_values,
    seed = 42
  )
}

# Chain 1: Unrecoverable error evaluating the log probability at the initial
# value. Chain 1: mismatch in number dimensions declared and found in context;
# processing stage=parameter initialization; variable name=Gamma; dims
# declared=(1,21); dims found=(4000,1,21) Error in sampler$call_sampler(c(args,
# dotlist)) : mismatch in number dimensions declared and found in context;
# processing stage=parameter initialization; variable name=Gamma; dims
# declared=(1,21); dims found=(4000,1,21)

# Check that fixing values is working (for full posterior).
k <- 2
beta_ids <- ensemble_fit[[k]] %>%
  gather_draws(Beta[r, i]) %>% 
  filter(r <= 5) %>% 
  unite(.variable, .variable, r, i) %>%
  distinct(.variable) %>%
  mutate(id = row_number()) %>%
  select(.variable, id)

ensemble_fit[[k]] %>%
  gather_draws(Beta[r, i]) %>%
  unite(.variable, .variable, r, i) %>%
  right_join(beta_ids) %>%
  ggplot(aes(x = .value, y = .variable)) +
  stat_halfeye(.width = .95) +
  facet_wrap(
    ~ .variable,
    ncol = dim(data$X)[4],
    scales = "free"
  )

ggsave(
  "marginals_check.png",
  path = here::here("Figures"),
  width = 30, height = 10, units = "in"
)

which(ind_ana[k,]==1)

# Save data and ensemble output.
data$ensemble_fit <- ensemble_fit
write_rds(
  data,
  here::here("Output", "fit_pathology-none.rds")
  # here::here("Output", "fit_vb_pathology-none.rds")
)
