# Preamble ----------------------------------------------------------------
# Load packages.
library(tidyverse)
library(cmdstanr)
library(posterior)
# library(bayesplot)
library(tidybayes)

# Set the simulation seed.
set.seed(42)

# Load data and draw a sample of ensembles of size nmember.
if (ind_sim == 1) data <- read_rds(here::here("data", str_c("sim_", file_id, ".rds")))
if (ind_emp == 1) data <- read_rds(here::here("data", str_c("emp_", data_id, ".rds")))
data$train_Z <- matrix(rep(1, nrow(data$train_Y)), ncol = 1)

######################
# The following mat_* code doesn't work for heterogeneous pathologies.
# Confirm all changes continue to work for homogeneous pathologies.
# Temporary fix:
if (ind_test == 1) {
  # If it's a test, we can't draw from mat_ana and mat_screen -- we use the whole thing.
  mat_ana <- data$mat_ana
  mat_screen <- data$mat_screen
}

# # These need to handle both matrices (homogeneous pathologies across members) and arrays
# # (heterogeneous pathologies across members).
# mat_ana <- data$mat_ana[sample(nrow(data$mat_ana), nmember),]
# mat_screen <- data$mat_screen[sample(nrow(data$mat_screen), nmember),]
# # mat_resp <- data$mat_resp[sample(nrow(data$mat_resp), nmember),]
######################

# Run HMNL to Initialize Ensemble -----------------------------------------
# Pathfinder alone could be used to initialize, if needed.
if (!file.exists(here::here("output", str_c("hmnl-fit_", ifelse(ind_sim == 1, file_id, data_id), ".rds")))) {
  # Specify data as a list.
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
  
  # Compile and estimate the model.
  hmnl <- cmdstan_model(here::here("code", "src", "hmnl.stan"))
  hmnl_fit <- hmnl$sample(
    data = stan_data,
    seed = 42,
    chains = 4,
    parallel_chains = 4
  )

  # Save HMNL fit.
  # write_rds(hmnl_fit, here::here("output", str_c("hmnl-fit_", file_id, ".rds")))
  hmnl_fit$save_object(here::here("output", str_c("hmnl-fit_", file_id, ".rds")))
} else {
  # Read HMNL fit.
  if (ind_sim == 1) hmnl_fit <- read_rds(here::here("output", str_c("hmnl-fit_", file_id, ".rds")))
  if (ind_emp == 1) hmnl_fit <- read_rds(here::here("output", str_c("hmnl-fit_", data_id, ".rds")))
}

# Use posteriors to construct priors for the ensemble.
hmnl_draws <- hmnl_fit$draws(format = "df")

Gamma_mean <- hmnl_draws |> 
  select(contains("Gamma")) |> 
  pivot_longer(everything()) |> 
  summarize(mean = mean(value)) |> 
  pull()
  
Gamma_scale <- hmnl_draws |> 
  select(contains("Gamma")) |> 
  pivot_longer(everything()) |> 
  summarize(sd = sd(value)) |> 
  pull()

Omega_shape <- hmnl_draws |> 
  select(contains("Omega")) |> 
  pivot_longer(everything()) |> 
  summarize(mean = mean(value)) |> 
  pull()

tau_mean <- hmnl_draws |> 
  select(contains("tau")) |> 
  pivot_longer(everything()) |> 
  summarize(mean = mean(value)) |> 
  pull()

tau_scale <- hmnl_draws |> 
  select(contains("tau")) |> 
  pivot_longer(everything()) |> 
  summarize(sd = sd(value)) |> 
  pull()

# Run HMNL Ensemble -------------------------------------------------------
# Pathfinder could be run in place of the full posterior, perhaps without any initialization.
# Specify data as a list of lists.
stan_data_list <- vector(mode = "list", length = nmember)
for (k in 1:nmember) {
  # Temporary full dataset.
  train_Y1 <- data$train_Y
  train_X1 <- data$train_X
  train_Z1 <- data$train_Z
  
  # # Reconstruct data for each ensemble member for respondent quality.
  # if (ind_resp == 1) {
  #   # Pull respondents for the specific member.
  #   mat_vec <- mat_resp[k,]
  #   for (i in 1:length(mat_vec)) {
  #     train_Y1[i,] <- data$train_Y[mat_vec[i],]
  #     train_X1[i,,,] <- data$train_X[mat_vec[i],,,]
  #     train_Z1[i,] <- data$train_Z[mat_vec[i],]
  #   }
  #   
  #   # # Overwrite the original data. WHY?! DON'T DO THIS.
  #   # data$train_Y <- train_Y1
  #   # data$train_X <- train_X1
  #   # data$train_Z <- train_Z1
  # }
  
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
    
    ##############
    # These need to handle both matrices (homogeneous pathologies across members) and arrays
    # (heterogeneous pathologies across members).
    # NO -- ARRAYS FOR BOTH.
    ##############
    
    mat_ana = mat_ana,         # Matrix of ensemble indicators for ANA.
    mat_screen = mat_screen    # Matrix of ensemble indicators for screening.
    
    ##############
    # Need to get matrices to work for homogeneous and heterogeneous pathologies.
    ##############
    
    # mat_ana = matrix(mat_ana, nrow = 1),         # Matrix of ensemble indicators for ANA.
    # mat_screen = matrix(mat_screen, nrow = 1)    # Matrix of ensemble indicators for screening.
  )
  
  stan_data_list[[k]] <- stan_data
}

fit_extract_average <- function(stan_data) {
  # # Estimate with VB.
  # fit <- rstan::vb(
  #   stan_model(here::here("Code", "src", "hmnl_ensemble.stan")),
  #   data = stan_data,
  #   init = 0,
  #   seed = 42,
  #   output_samples = 100
  # )
  
  # Compile and estimate the model.
  hmnl_ensemble <- cmdstan_model(here::here("code", "src", "hmnl_ensemble_01.stan"))
  # hmnl_ensemble <- cmdstan_model(here::here("code", "src", "hmnl_ensemble_02.stan"))
  fit <- hmnl_ensemble$sample(
    data = stan_data,
    seed = 42,
    chains = 1,
    thin = 10
    # chains = 4,
    # parallel_chains = 4
  )

  # Extract the posterior draws for Gamma, Sigma, and log_lik.
  # draws <- fit$draws(format = "df", variables = c("Gamma", "Sigma", "log_lik"))
  draws <- fit$draws(format = "df", variables = c("Beta", "Gamma", "Sigma", "log_lik"))
  
  # Compute posterior means.
  ensemble_draws <- NULL
  ensemble_draws$Beta <- draws |> subset_draws(variable = "Beta") #|> summarize_draws("mean")
  ensemble_draws$Gamma <- draws |> subset_draws(variable = "Gamma") |> summarize_draws("mean")
  ensemble_draws$Sigma <- draws |> subset_draws(variable = "Sigma") |> summarize_draws("mean")
  ensemble_draws$log_lik <- draws |> subset_draws(variable = "log_lik")
  
  return(ensemble_draws)
}

ensemble_draws <- parallel::mclapply(stan_data_list, fit_extract_average, mc.cores = parallel::detectCores())

# Save ensemble fit.
ensemble_fit <- list(mat_ana = mat_ana, mat_screen = mat_screen, ensemble_draws = ensemble_draws)
write_rds(ensemble_fit, here::here("output", str_c("ensemble-fit_", file_id, "_", nmember, ".rds")))

# rempBayes(lgtdata = stan_data) # WHAT IS THIS?

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
#   #   stan_model(here::here("Code", "src", "hmnl_ensemble.stan")),
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
#     here::here("Code", "src", "hmnl_ensemble.stan"),
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
