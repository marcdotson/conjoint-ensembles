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
ind_none <- 0        # Indicates no pathologies.
ind_ana <- 1        # Indicates attribute non-attendance.
ind_screen <- 0     # Indicates screening.
ind_ana_screen <- 0 # Indicates attribute non-attendance and screening.
ind_real <- 0       # Indicates ____ data.

if (ind_non == 1) file_name <- "none"
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

initial_fit <- stan(
  here::here("Code", "Source", "hmnl.stan"),
  data = stan_data,
  seed = 42
)

# Save initial fit output.
write_rds(
  initial_fit,
  here::here("Output", str_c("initial-fit_", file_name, ".rds"))
)

# Load initial fit output.
initial_fit <- read_rds(here::here("Output", str_c("initial-fit_", file_name, ".rds")))

# Construct initial values.
initial_draws <- extract(initial_fit, pars = c("Gamma", "Omega", "tau", "Delta", "Beta"))
initial_values <- vector(mode = "list", length = 1)

initial_values[[1]]$Gamma <- initial_draws$Gamma %>% 
  array_branch(margin = c(2, 3)) %>% 
  map(mean) %>% 
  matrix(nrow = dim(initial_draws$Gamma)[2], ncol = dim(initial_draws$Gamma)[3])

initial_values[[1]]$Omega <- initial_draws$Omega %>% 
  array_branch(margin = c(2, 3)) %>% 
  map(mean) %>% 
  matrix(nrow = dim(initial_draws$Omega)[2], ncol = dim(initial_draws$Omega)[3])

initial_values[[1]]$tau <- initial_draws$tau %>% 
  array_branch(margin = 2) %>% 
  map(mean) %>% 
  matrix(nrow = dim(initial_draws$Gamma)[2])

initial_values[[1]]$Delta <- initial_draws$Delta %>% 
  array_branch(margin = c(2, 3)) %>% 
  map(mean) %>% 
  matrix(nrow = dim(initial_draws$Delta)[2], ncol = dim(initial_draws$Delta)[3])

initial_values[[1]]$Beta <- initial_draws$Beta %>% 
  array_branch(margin = c(2, 3)) %>% 
  map(mean) %>% 
  matrix(nrow = dim(initial_draws$Beta)[2], ncol = dim(initial_draws$Beta)[3])

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
  
  ensemble_fit[[k]] <- vb(
    stan_model(here::here("Code", "Source", "hmnl_ensemble.stan")),
    data = stan_data,
    init = 0,
    # init = initial_values,
    seed = 42
  )
  
  # ensemble_fit[[k]] <- stan(
  #   here::here("Code", "Source", "hmnl_ensemble.stan"),
  #   data = stan_data,
  #   init = initial_values,
  #   seed = 42
  # )
}

# Save data and ensemble output.
data$ensemble_fit <- ensemble_fit
write_rds(
  data,
  here::here("Output", str_c("fit-vb_", file_name, ".rds"))
)

# Load data and ensemble output.
data <- read_rds(here::here("Output", str_c("fit-vb_", file_name, ".rds")))

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

