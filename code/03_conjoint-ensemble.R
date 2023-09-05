# Preamble ----------------------------------------------------------------
# Load packages.
library(tidyverse)
# library(rstan) # REPLACE WITH CMDSTANR
library(cmdstanr)
library(posterior)
library(bayesplot)
library(tidybayes)

# # Set Stan and future options.
# options(mc.cores = parallel::detectCores())
# rstan_options(auto_write = FALSE)

# Set the simulation seed.
set.seed(42)

# Load data and draw a sample of ensembles of size nmember.
if (ind_sim == 1) data <- read_rds(here::here("data", str_c("sim_", file_id, ".rds")))
if (ind_emp == 1) data <- read_rds(here::here("data", str_c("emp_", data_id, ".rds")))
data$train_Z <- matrix(rep(1, nrow(data$train_Y)), ncol = 1)
mat_ana <- data$mat_ana[sample(nrow(data$mat_ana), nmember),]
mat_screen <- data$mat_screen[sample(nrow(data$mat_screen), nmember),]
mat_resp <- data$mat_resp[sample(nrow(data$mat_resp), nmember),]

# Run HMNL to Initialize Ensemble -----------------------------------------
#################################
# Pathfinder alone should be used to initalize (if even needed).
#################################
if (!file.exists(here::here("output", str_c("hmnl-fit_", ifelse(ind_sim == 1, file_id, data_id), ".rds")))) {
  stan_data <- list(
    R = dim(data$train_X)[1], # Number of respondents.
    S = dim(data$train_X)[2], # Number of choice tasks.
    A = dim(data$train_X)[3], # Number of choice alternatives.
    I = dim(data$train_X)[4], # Number of observation-level covariates.
    J = ncol(data$train_Z),   # Number of population-level covariates.
    
    Gamma_mean = 0,           # Mean of population-level means.
    Gamma_scale = 5,          # Scale of population-level means.
    Omega_shape = 2,          # Shape of population-level scale.
    tau_mean = 0,             # Mean of population-level scale.
    tau_scale = 5,            # Scale of population-level scale.
    
    Y = data$train_Y,         # Matrix of observations.
    X = data$train_X,         # Array of observation-level covariates.
    Z = data$train_Z          # Matrix of population-level covariates.
  )
  
  # hmnl_fit <- stan(
  #   here::here("code", "source", "hmnl.stan"),
  #   data = stan_data,
  #   seed = 42
  # )
  
  # Save HMNL fit.
  write_rds(hmnl_fit, here::here("output", str_c("hmnl-fit_", file_id, ".rds")))
} else {
  # Read HMNL fit.
  if (ind_sim == 1) hmnl_fit <- read_rds(here::here("output", str_c("hmnl-fit_", file_id, ".rds")))
  if (ind_emp == 1) hmnl_fit <- read_rds(here::here("output", str_c("hmnl-fit_", data_id, ".rds")))
}

# # Use posteriors to construct priors.
# hmnl_draws <- extract(hmnl_fit, pars = c("Gamma", "Omega", "tau", "Delta"))
# Gamma_mean <- mean(hmnl_draws$Gamma)
# Gamma_scale <- sqrt(var(hmnl_draws$Gamma))
# Omega_shape <- mean(hmnl_draws$Omega)
# tau_mean <- mean(hmnl_draws$tau)
# tau_scale <- sqrt(var(as.vector(hmnl_draws$tau)))

# Run HMNL Ensemble -------------------------------------------------------
# # NOT PARALLELIZED
# temp <- vector(mode = "list", length = nmember)
# for (k in 1:nmember) {
#   # Temporary full dataset.
#   train_Y1 = data$train_Y
#   train_X1 = data$train_X
#   train_Z1 <- data$train_Z
#   
#   # Reconstruct data for each ensemble member for respondent quality.
#   if (ind_resp == 1) {
#     # Pull respondents for the specific member.
#     mat_vec <- mat_resp[k,]
#     for (i in 1:length(mat_vec)) {
#       train_Y1[i,] <- data$train_Y[mat_vec[i],]
#       train_X1[i,,,] <- data$train_X[mat_vec[i],,,]
#       train_Z1[i,] <- data$train_Z[mat_vec[i],]
#     }
# 
#     # # Overwrite the original data. WHY?! DON'T DO THIS.
#     # data$train_Y <- train_Y1
#     # data$train_X <- train_X1
#     # data$train_Z <- train_Z1
#   }
# 
#   stan_data <- list(
#     R = dim(data$train_X)[1],  # Number of respondents.
#     S = dim(data$train_X)[2],  # Number of choice tasks.
#     A = dim(data$train_X)[3],  # Number of choice alternatives.
#     I = dim(data$train_X)[4],  # Number of observation-level covariates.
#     J = ncol(data$train_Z),    # Number of population-level covariates.
#     K = nmember,               # Number of members of the ensemble.
#     k = k,                     # Ensemble member number.
# 
#     Gamma_mean = Gamma_mean,   # Mean of population-level means.
#     Gamma_scale = Gamma_scale, # Scale of population-level means.
#     Omega_shape = 2,           # Shape of population-level scale.
#     tau_mean = tau_mean,       # Mean of population-level scale.
#     tau_scale = tau_scale,     # Scale of population-level scale.
# 
#     # Y = data$train_Y,          # Matrix of observations.
#     # X = data$train_X,          # Array of observation-level covariates.
#     # Z = data$train_Z,          # Matrix of population-level covariates.
#     Y = train_Y1,          # Matrix of observations.
#     X = train_X1,          # Array of observation-level covariates.
#     Z = train_Z1,          # Matrix of population-level covariates.
# 
#     ind_ana = ind_ana,         # Flag indicating attribute non-attendance.
#     ind_screen = ind_screen,   # Flag indicating screening.
#     mat_ana = mat_ana,         # Matrix of ensemble indicators for ANA.
#     mat_screen = mat_screen    # Matrix of ensemble indicators for screening.
#   )
# 
#   # init_fun <- function(...) {
#   #   list(
#   #     Gamma = apply(hmnl_draws$Gamma, c(2, 3), mean),
#   #     Omega = apply(hmnl_draws$Omega, c(2, 3), mean),
#   #     tau = apply(hmnl_draws$tau, 2, mean),
#   #     Delta = apply(hmnl_draws$Delta, c(2, 3), mean)
#   #   )
#   #   # matrix[J, I] Gamma;                // Matrix of population-level hyperparameters.
#   #   # corr_matrix[I] Omega;              // Population model correlation matrix hyperparameters.
#   #   # vector<lower = 0>[I] tau;          // Population model vector of scale hyperparameters.
#   #   # matrix[R, I] Delta;                // Matrix of non-centered observation-level parameters.
#   # }
# 
#   # # Estimate with VB.
#   # fit <- rstan::vb(
#   #   stan_model(here::here("Code", "Source", "hmnl_ensemble.stan")),
#   #   data = stan_data,
#   #   init = 0,
#   #   # init = init_fun,
#   #   seed = 42,
#   #   # eval_elbo = 200,
#   #   # elbo_samples = 200,
#   #   output_samples = 100
#   # )
#   
#   # Estimate with full posterior sampling.
#   fit <- stan(
#     here::here("Code", "Source", "hmnl_ensemble.stan"),
#     data = stan_data,
#     chains = 1,
#     thin = 10,
#     seed = 42
#   )
# 
#   # Extract the posterior draws for Gamma, Sigma, and log_lik.
#   draws <- rstan::extract(fit, pars = c("Gamma", "Sigma", "log_lik"))
# 
#   # Compute posterior means.
#   ensemble_draws <- NULL
#   ensemble_draws$Gamma <- apply(draws$Gamma, c(2, 3), mean)
#   ensemble_draws$Sigma <- apply(draws$Sigma, c(2, 3), mean)
#   ensemble_draws$log_lik <- draws$log_lik
# 
#   temp[[k]] <- ensemble_draws
# }



# PARALLELIZED
# Compute the conjoint ensemble.
stan_data_list <- vector(mode = "list", length = nmember)
for (k in 1:nmember) {
  # Temporary full dataset.
  train_Y1 = data$train_Y
  train_X1 = data$train_X
  train_Z1 <- data$train_Z
  
  # Reconstruct data for each ensemble member for respondent quality.
  if (ind_resp == 1) {
    # Pull respondents for the specific member.
    mat_vec <- mat_resp[k,]
    for (i in 1:length(mat_vec)) {
      train_Y1[i,] <- data$train_Y[mat_vec[i],]
      train_X1[i,,,] <- data$train_X[mat_vec[i],,,]
      train_Z1[i,] <- data$train_Z[mat_vec[i],]
    }
    
    # # Overwrite the original data. WHY?! DON'T DO THIS.
    # data$train_Y <- train_Y1
    # data$train_X <- train_X1
    # data$train_Z <- train_Z1
  }
  
  stan_data <- list(
    R = dim(data$train_X)[1],  # Number of respondents.
    S = dim(data$train_X)[2],  # Number of choice tasks.
    A = dim(data$train_X)[3],  # Number of choice alternatives.
    I = dim(data$train_X)[4],  # Number of observation-level covariates.
    J = ncol(data$train_Z),    # Number of population-level covariates.
    K = nmember,               # Number of members of the ensemble.
    k = k,                     # Ensemble member number.
    
    Gamma_mean = Gamma_mean,   # Mean of population-level means.
    Gamma_scale = Gamma_scale, # Scale of population-level means.
    Omega_shape = 2,           # Shape of population-level scale.
    tau_mean = tau_mean,       # Mean of population-level scale.
    tau_scale = tau_scale,     # Scale of population-level scale.
    
    # Y = data$train_Y,          # Matrix of observations.
    # X = data$train_X,          # Array of observation-level covariates.
    # Z = data$train_Z,          # Matrix of population-level covariates.
    Y = train_Y1,          # Matrix of observations.
    X = train_X1,          # Array of observation-level covariates.
    Z = train_Z1,          # Matrix of population-level covariates.
    
    ind_ana = ind_ana,         # Flag indicating attribute non-attendance.
    ind_screen = ind_screen,   # Flag indicating screening.
    mat_ana = mat_ana,         # Matrix of ensemble indicators for ANA.
    mat_screen = mat_screen    # Matrix of ensemble indicators for screening.
  )
  
  stan_data_list[[k]] <- stan_data
}

fit_extract_average <- function(stan_data) {
  # # Estimate with VB.
  # fit <- rstan::vb(
  #   stan_model(here::here("Code", "Source", "hmnl_ensemble.stan")),
  #   data = stan_data,
  #   init = 0,
  #   seed = 42,
  #   output_samples = 100
  # )
  
  # # Estimate with full posterior sampling.
  # fit <- stan(
  #   here::here("code", "source", "hmnl_ensemble.stan"),
  #   data = stan_data,
  #   chains = 1,
  #   thin = 10,
  #   seed = 42
  # )
  
  # # Extract the posterior draws for Gamma, Sigma, and log_lik.
  # draws <- rstan::extract(fit, pars = c("Gamma", "Sigma", "log_lik"))
  
  # Compute posterior means.
  ensemble_draws <- NULL
  ensemble_draws$Gamma <- apply(draws$Gamma, c(2, 3), mean)
  ensemble_draws$Sigma <- apply(draws$Sigma, c(2, 3), mean)
  ensemble_draws$log_lik <- draws$log_lik
  
  return(ensemble_draws)
}

# ensemble_draws <- parallel::mclapply(stan_data_list, fit_extract_average, mc.cores = parallel::detectCores())

# Save ensemble fit.
ensemble_fit <- list(mat_ana = mat_ana, mat_screen = mat_screen, ensemble_draws = ensemble_draws)
write_rds(ensemble_fit, here::here("output", str_c("ensemble-fit_", file_id, "_", nmember, ".rds")))

