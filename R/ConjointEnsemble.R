# R script to run a pathology treatment ensemble
# Preamble ----------------------------------------------------------------
# Load libraries.
library(loo)
library(rstan)

# Set Stan options.
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

nresp = 100
ntask = 10
nalts = 4
nlvls = 12
ncovs = 1

# Save data set variables in a list.
dataset_vars <- list(R = nresp, T = ntask, A = nalts, L = nlvls, C = ncovs)
stan_data <- stan("./STAN/generate_data.stan",
                  data = dataset_vars,
                  iter = 1,
                  chains = 1,
                  warmup = 0,
                  refresh = 1,
                  seed = 1750532,
                  algorithm = 'Fixed_param')

# Extract stan data from fit
simX <- extract(stan_data, permuted = TRUE)$X
simY <- extract(stan_data, permuted = TRUE)$Y
simZ <- extract(stan_data, permuted = TRUE)$Z

# store in a list to pass to the sampling method
dataset_vars <- list(R = nresp,
                     T = ntask,
                     A = nalts,
                     L = nlvls,
                     C = ncovs,
                     X = as.array(simX[1,,,,]),
                     Y = as.array(simY[1,,]),
                     Z = as.matrix(simZ[1,,]))

# Use the following Stan models
model_list <- c("./STAN/mnl_vanilla.stan",
                "./STAN/mnl_vanilla.stan",
                "./STAN/mnl_vanilla.stan",
                "./STAN/mnl_vanilla.stan",
                "./STAN/mnl_fhorseshoe.stan",
                "./STAN/mnl_fhorseshoe.stan",
                "./STAN/mnl_fhorseshoe.stan",
                "./STAN/mnl_fhorseshoe.stan",
                "./STAN/mnl_leakyrelu.stan",
                "./STAN/mnl_leakyrelu.stan",
                "./STAN/mnl_leakyrelu.stan",
                "./STAN/mnl_leakyrelu.stan")

# Run the ensemble via loo package
fit_list <- list()
K <- length(model_list)
for (k in 1:K){
    # Fit the k-th model with Stan
    fit <- stan(model_list[[k]],
                data = dataset_vars,
                iter=500,
                chains = 2,
                control = list(adapt_delta = .9, max_treedepth = 3))
    fit_list[[k]] <- fit
}

log_lik_list <- lapply(fit_list, extract_log_lik)

r_eff_list <- lapply(fit_list, function(x) {
  ll_array <- extract_log_lik(x, merge_chains = FALSE)
  relative_eff(exp(ll_array))
})

W <- loo_model_weights(log_lik_list,
                       method = 'stacking',
                       r_eff_list = r_eff_list,
                       optim_control = list(reltol = 1e-10))
print(W)
