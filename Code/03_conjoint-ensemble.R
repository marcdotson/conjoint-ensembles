# Preamble ----------------------------------------------------------------
# Load packages.
library(tidyverse)
library(rstan)
# library(bayesplot)
# library(tidybayes)
# library(loo)

# Set Stan options.
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Load Data ---------------------------------------------------------------
ind_none <- 1       # Indicates no pathologies.
ind_ana <- 0        # Indicates attribute non-attendance.
ind_screen <- 0     # Indicates screening.
ind_ana_screen <- 0 # Indicates attribute non-attendance and screening.
ind_real <- 0       # Indicates ____ data.

if (ind_none == 1) file_name <- "none"
if (ind_ana == 1) file_name <- "ana"
if (ind_screen == 1) file_name <- "screen"
if (ind_ana_screen == 1) file_name <- "ana-screen"
if (ind_real == 1) file_name <- "design"

data <- read_rds(here::here("Data", str_c("sim_", file_name, ".rds")))
Y <- data$train_Y
X <- data$train_X
Z <- matrix(rep(1, nrow(Y)), ncol = 1)
mat_ana <- data$mat_ana
mat_screen <- data$mat_screen

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

hmnl_fit <- stan(
  here::here("Code", "Source", "hmnl.stan"),
  data = stan_data,
  seed = 42
)

# Save hmnl-fit output.
write_rds(
  hmnl_fit,
  here::here("Output", str_c("hmnl-fit_", file_name, ".rds"))
)

# Load hmnl-fit output.
hmnl_fit <- read_rds(here::here("Output", str_c("hmnl-fit_", file_name, ".rds")))

# Construct initial values.
# parameters <- c("Gamma", "Omega", "tau", "Delta", "Beta")
parameters <- c("Gamma", "Omega", "tau", "Delta")
hmnl_draws <- extract(hmnl_fit, pars = parameters)
initial_values <- vector(mode = "list", length = 1)
for (i in 1:length(parameters)) {
  parameter <- parameters[i]
  if (parameter != "tau") {
    initial_values[[1]][[parameter]] <- hmnl_draws[parameter][[1]] %>% 
      array_branch(margin = c(2, 3)) %>% 
      map(mean) %>% 
      unlist() %>% 
      matrix(nrow = dim(hmnl_draws[parameter][[1]])[2], ncol = dim(hmnl_draws[parameter][[1]])[3])
  }
  if (parameter == "tau") {
    initial_values[[1]][[parameter]] <- hmnl_draws[parameter][[1]] %>% 
      array_branch(margin = 2) %>% 
      map(mean) %>% 
      unlist() %>% 
      as.vector()
  }
}
for (j in 2:4) initial_values[[j]] <- initial_values[[1]]

# Ensemble Calibration ----------------------------------------------------
K <- nrow(mat_ana)
ensemble_fit <- vector(mode = "list", length = K)
for (k in 1:K) {
  stan_data <- list(
    R = dim(X)[1],          # Number of respondents.
    S = dim(X)[2],          # Number of choice tasks.
    A = dim(X)[3],          # Number of choice alternatives.
    I = dim(X)[4],          # Number of observation-level covariates.
    J = ncol(Z),            # Number of population-level covariates.
    K = K,                  # Number of members of the ensemble.
    k = k,                  # Ensemble member number.
    
    Gamma_mean = 0,         # Mean of population-level means.
    Gamma_scale = 5,        # Scale of population-level means.
    Omega_shape = 2,        # Shape of population-level scale.
    tau_mean = 0,           # Mean of population-level scale.
    tau_scale = 5,          # Scale of population-level scale.
    
    Y = Y,                  # Matrix of observations.
    X = X,                  # Array of observation-level covariates.
    Z = Z,                  # Matrix of population-level covariates.
    mat_ana = mat_ana       # Matrix of ensemble indicators for ANA.
    # mat_screen = mat_screen # Matrix of ensemble indicators for screening.
  )
  
  # Set initial values by providing a list equal in length to the number of chains.
  # The elements of this list should themselves be named lists, where each of
  # these named lists has the name of a parameter and is used to specify the
  # initial values for that parameter for the corresponding chain.
  
  # ensemble_fit[[k]] <- vb(
  #   stan_model(here::here("Code", "Source", "hmnl_ensemble.stan")),
  #   data = stan_data,
  #   # init = 0,
  #   init = initial_values,
  #   seed = 42
  # )
  
  ensemble_fit[[k]] <- stan(
    here::here("Code", "Source", "hmnl_ensemble.stan"),
    data = stan_data,
    init = initial_values,
    seed = 42
  )
}

# Save data and ensemble output.
data$ensemble_fit <- ensemble_fit
write_rds(
  data,
  # here::here("Output", str_c("ensemble-fit_vb_", file_name, ".rds"))
  here::here("Output", str_c("ensemble-fit_", file_name, ".rds"))
)

# Load data and ensemble output.
data <- read_rds(here::here("Output", str_c("ensemble-fit_vb_", file_name, ".rds")))

# Post-hoc extraction and reassembly of ensemble output.
for (k in 1:nrow(mat_ana)) {
  data$ensemble_fit[[k]] <- extract(data$ensemble_fit[[k]], pars = c("log_lik", "Beta"))
}
write_rds(
  data,
  here::here("Output", str_c("reduced_ensemble-fit_vb_", file_name, ".rds"))
)

# # Check that fixing values is working (for full posterior).
# k <- 2
# beta_ids <- ensemble_fit[[k]] %>%
#   gather_draws(Beta[r, i]) %>% 
#   filter(r <= 5) %>% 
#   unite(.variable, .variable, r, i) %>%
#   distinct(.variable) %>%
#   mutate(id = row_number()) %>%
#   select(.variable, id)
# 
# ensemble_fit[[k]] %>%
#   gather_draws(Beta[r, i]) %>%
#   unite(.variable, .variable, r, i) %>%
#   right_join(beta_ids) %>%
#   ggplot(aes(x = .value, y = .variable)) +
#   stat_halfeye(.width = .95) +
#   facet_wrap(
#     ~ .variable,
#     ncol = dim(data$X)[4],
#     scales = "free"
#   )
# 
# ggsave(
#   "marginals_check.png",
#   path = here::here("Figures"),
#   width = 30, height = 10, units = "in"
# )
# 
# which(ind_ana[k,]==1)

