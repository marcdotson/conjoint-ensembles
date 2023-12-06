# Run the conjoint ensemble using clever randomization.
data <- read_rds(here::here("data", str_c(data_id, "_", patho_id, ".rds")))

# Initialize the ensemble by estimating a full HMNL and then using the posteriors of the full
# HMNL as the hyperprior values for the complete ensemble. Use the full HMNL to compare post-hoc
# modification of the betas rather than a complete ensemble.
if (!file.exists(here::here("output", str_c("hmnl-fit_", data_id, "_", patho_id, ".rds")))) {
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
  hmnl_fit$save_object(here::here("output", str_c("hmnl-fit_", data_id, "_", patho_id, ".rds")))
} else {
  # Read HMNL fit.
  hmnl_fit <- read_rds(here::here("output", str_c("hmnl-fit_", data_id, "_", patho_id, ".rds")))
}

# Extract posteriors draws and construct hyperpriors for the ensemble.
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

# Specify data as a list of lists.
stan_data_list <- vector(mode = "list", length = nmember)
for (member in 1:nmember) {
  stan_data <- list(
    R = dim(data$train_X)[1],                   # Number of respondents.
    S = dim(data$train_X)[2],                   # Number of choice tasks.
    A = dim(data$train_X)[3],                   # Number of choice alternatives.
    I = dim(data$train_X)[4],                   # Number of observation-level covariates.
    J = ncol(data$train_Z),                     # Number of population-level covariates.
              
    Gamma_mean = Gamma_mean,                    # Mean of population-level means.
    Gamma_scale = Gamma_scale,                  # Scale of population-level means.
    Omega_shape = Omega_shape,                  # Shape of population-level scale.
    tau_mean = tau_mean,                        # Mean of population-level scale.
    tau_scale = tau_scale,                      # Scale of population-level scale.
              
    Y = data$train_Y,                           # Array of observations.
    X = data$train_X,                           # Array of observation-level covariates.
    Z = data$train_Z,                           # Matrix of population-level covariates.
    array_ana = data$array_ana[,,member],       # Array of ensemble indicators for ANA.
    array_screen = data$array_screen[,,member], # Array of ensemble indicators for screening.
    
    ##########################
    # C'mon dimensions...
    ##########################
    
    array_qual = as.matrix(data$array_qual[,,member])      # Array of ensemble indicators for respondent quality.
  )
  
  stan_data_list[[member]] <- stan_data
}

##########################
# Do I need to compile every time?
# Do I need to pass a vector of arguments to parallel::mclapply()?
##########################

# Specify a function to fit each ensemble, extract the posterior draws, and average them.
fit_extract_average <- function(stan_data) {
  # Estimate the model.
  hmnl_ensemble <- cmdstan_model(here::here("code", "src", "hmnl_ensemble.stan"))
  fit <- hmnl_ensemble$sample(
    data = stan_data,
    seed = 42,
    chains = 1,
    thin = 10
  )

  # Extract the posterior draws for Gamma, Sigma, and log_lik.
  draws <- fit$draws(format = "df", variables = c("Beta", "Gamma", "Sigma", "log_lik"))
  
  # Compute posterior means.
  ensemble_draws <- NULL
  ensemble_draws$Beta <- draws |> subset_draws(variable = "Beta") #|> summarize_draws("mean")
  ensemble_draws$Gamma <- draws |> subset_draws(variable = "Gamma") |> summarize_draws("mean")
  ensemble_draws$Sigma <- draws |> subset_draws(variable = "Sigma") |> summarize_draws("mean")
  ensemble_draws$log_lik <- draws |> subset_draws(variable = "log_lik")
  
  return(ensemble_draws)
}

# Fit the ensemble.
ensemble_draws <- parallel::mclapply(stan_data_list, fit_extract_average, mc.cores = parallel::detectCores())

# Save the ensemble fit.
write_rds(ensemble_draws, here::here("output", str_c("ensemble-fit_", data_id, "_", patho_id, "_", nmember, ".rds")))

