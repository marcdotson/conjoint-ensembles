# Preamble ----------------------------------------------------------------
# Load packages.
library(tidyverse)
library(rstan)

# Set Stan options.
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Load Data ---------------------------------------------------------------
ind_none <- 0       # Indicates no pathologies.
ind_ana <- 1        # Indicates attribute non-attendance.
ind_screen <- 0     # Indicates screening.
ind_ana_screen <- 0 # Indicates attribute non-attendance and screening.

if (ind_none == 1) file_name <- "none"
if (ind_ana == 1) file_name <- "ana"
if (ind_screen == 1) file_name <- "screen"
if (ind_ana_screen == 1) file_name <- "ana-screen"

# data <- read_rds(here::here("Data", str_c("sim_", file_name, ".rds")))
data <- read_rds(here::here("Data", str_c("sim_", file_name, "_01.rds")))
Y <- data$train_Y
X <- data$train_X
Z <- matrix(rep(1, nrow(Y)), ncol = 1)
mat_ana <- data$mat_ana
mat_screen <- data$mat_screen

# Run HMNL to Initialize Ensemble -----------------------------------------
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

hmnl_fit <- stan(
  here::here("Code", "Source", "hmnl.stan"),
  data = stan_data,
  seed = 42
)

# Save HMNL fit.
write_rds(hmnl_fit, here::here("Output", str_c("hmnl-fit_", file_name, ".rds")))

# Load HMNL fit.
# hmnl_fit <- read_rds(here::here("Output", str_c("hmnl-fit_", file_name, ".rds")))
hmnl_fit <- read_rds(here::here("Output", str_c("hmnl-fit_", file_name, "_01.rds")))

# Use posteriors to construct priors.
hmnl_draws <- extract(hmnl_fit, pars = c("Gamma", "Omega", "tau"))
Gamma_mean <- mean(hmnl_draws$Gamma)
Gamma_scale <- sqrt(var(hmnl_draws$Gamma))
Omega_shape <- mean(hmnl_draws$Omega)
tau_mean <- mean(hmnl_draws$tau)
tau_scale <- sqrt(var(as.vector(hmnl_draws$tau)))
# parameters <- c("Gamma", "Omega", "tau", "Delta")
# hmnl_draws <- extract(hmnl_fit, pars = parameters)
# initial_values <- vector(mode = "list", length = 1)
# for (i in 1:length(parameters)) {
#   parameter <- parameters[i]
#   if (parameter != "tau") {
#     initial_values[[1]][[parameter]] <- hmnl_draws[parameter][[1]] %>% 
#       array_branch(margin = c(2, 3)) %>% 
#       map(mean) %>% 
#       unlist() %>% 
#       matrix(nrow = dim(hmnl_draws[parameter][[1]])[2], ncol = dim(hmnl_draws[parameter][[1]])[3])
#   }
#   if (parameter == "tau") {
#     initial_values[[1]][[parameter]] <- hmnl_draws[parameter][[1]] %>% 
#       array_branch(margin = 2) %>% 
#       map(mean) %>% 
#       unlist() %>% 
#       as.vector()
#   }
# }

# Run HMNL Ensemble -------------------------------------------------------
K <- nrow(mat_ana)
ensemble_fit <- vector(mode = "list", length = K)
for (k in 1:K) {
  stan_data <- list(
    R = dim(X)[1],             # Number of respondents.
    S = dim(X)[2],             # Number of choice tasks.
    A = dim(X)[3],             # Number of choice alternatives.
    I = dim(X)[4],             # Number of observation-level covariates.
    J = ncol(Z),               # Number of population-level covariates.
    K = K,                     # Number of members of the ensemble.
    k = k,                     # Ensemble member number.
    
    # Gamma_mean = 0,         # Mean of population-level means.
    # Gamma_scale = 5,        # Scale of population-level means.
    Omega_shape = 2,        # Shape of population-level scale.
    # tau_mean = 0,           # Mean of population-level scale.
    # tau_scale = 5,          # Scale of population-level scale.
    
    Gamma_mean = Gamma_mean,   # Mean of population-level means.
    Gamma_scale = Gamma_scale, # Scale of population-level means.
    # Omega_shape = Omega_shape, # Shape of population-level scale.
    tau_mean = tau_mean,       # Mean of population-level scale.
    tau_scale = tau_scale,     # Scale of population-level scale.
    
    Y = Y,                     # Matrix of observations.
    X = X,                     # Array of observation-level covariates.
    Z = Z,                     # Matrix of population-level covariates.
    mat_ana = mat_ana          # Matrix of ensemble indicators for ANA.
    # mat_screen = mat_screen    # Matrix of ensemble indicators for screening.
  )
  
  ensemble_fit[[k]] <- vb(
    stan_model(here::here("Code", "Source", "hmnl_ensemble.stan")),
    data = stan_data,
    init = 0,
    seed = 42
  )
  
  # ensemble_fit[[k]] <- stan(
  #   here::here("Code", "Source", "hmnl_ensemble.stan"),
  #   data = stan_data,
  #   init = initial_values,
  #   seed = 42
  # )
}

# Save ensemble fit.
# write_rds(ensemble_fit, here::here("Output", str_c("ensemble-fit_vb_", file_name, ".rds")))
write_rds(ensemble_fit, here::here("Output", str_c("ensemble-fit_vb_", file_name, "_01.rds")))

# Load ensemble fit.
ensemble_fit <- read_rds(here::here("Output", str_c("ensemble-fit_vb_", file_name, ".rds")))

