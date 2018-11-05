# R script to run a pathology treatment ensemble
# Preamble ----------------------------------------------------------------
# Load libraries.
Sys.setenv(USE_CXX14 = 1)
library(parallel)
library(loo)
library(rstan)

# Set Stan options.
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Save data set variables in a list.
dataset_vars <- list(R = 100, T = 10, A = 4, L = 12, C = 1)
stan_data <- stan("./STAN/generate_data.stan",
                  data = dataset_vars,
                  iter = 1,
                  chains = 1,
                  warmup = 0,
                  refresh = 1,
                  algorithm = 'Fixed_param')

# Extract stan data from fit
simX <- extract(stan_data, permuted = TRUE)$X
simY <- extract(stan_data, permuted = TRUE)$Y
simZ <- extract(stan_data, permuted = TRUE)$Z

# store in a list to pass to the sampling method
dataset_vars <- list(R = 100,
                     T = 10,
                     A = 4,
                     L = 12,
                     C = 1,
                     X = as.array(simX[1,,,,]),
                     Y = as.array(simY[1,,]),
                     Z = as.matrix(simZ[1,,]))

# Use the following Stan models
model_list <- c("./STAN/mnl_vanilla.stan",
                "./STAN/mnl_fhorseshoe.stan",
                "./STAN/mnl_leakyrelu.stan")

# Run the ensemble via loo package
log_lik_list <- list()
K <- 3
for (k in 1:K){
    # Fit the k-th model with Stan
    fit <- stan(model_list[[k]], data = dataset_vars, chains = 4)
    log_lik_list[[k]] <- extract_log_lik(fit)
}
model_weights <- loo_model_weights(log_lik_list, method="stacking")
print(model_weights)
