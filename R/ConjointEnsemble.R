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
niter = 500

# input data for generate_data.stan
train_vars <- list(R = nresp, T = ntask, A = nalts, L = nlvls, C = ncovs)

# Generate training data set
train_data <- stan("./STAN/generate_data.stan",
                   data = train_vars,
                   iter = 1,
                   chains = 1,
                   warmup = 0,
                   refresh = 1,
                   seed = 1750532,
                   algorithm = 'Fixed_param')
                  
# Extract training data from fit
simX <- extract(train_data, permuted = TRUE)$X
simY <- extract(train_data, permuted = TRUE)$Y
simZ <- extract(train_data, permuted = TRUE)$Z

# input data for the Stan base models
training_vars <- list(R = nresp,
                      T = ntask,
                      A = nalts,
                      L = nlvls,
                      C = ncovs,
                      X = as.array(simX[1,,,,]),
                      Y = as.array(simY[1,,]),
                      Z = as.matrix(simZ[1,,]))

# Use the following Stan models in the ensemble
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
#model_list <- c("./STAN/mnl_vanilla.stan",
#                "./STAN/mnl_fhorseshoe.stan",
#                "./STAN/mnl_leakyrelu.stan")


# Fit each base model to training data
fit_list <- list()
K <- length(model_list)
for (k in 1:K){
    # Fit the k-th model with Stan
    fit <- stan(model_list[[k]],
                data = training_vars,
                iter=niter,
                chains = 2,
                control = list(adapt_delta = .9, max_treedepth = 3))
    fit_list[[k]] <- fit
}


### COMPUTE STACKING WEIGHTS ###

# store each model's log likelihoods in a list
log_lik_list <- lapply(fit_list, extract_log_lik)

# relative effective sample size helps loo computation accuracy
r_eff_list <- lapply(fit_list,
                     function(x) {
                       ll_array <- extract_log_lik(x, merge_chains = FALSE)
                       relative_eff(exp(ll_array))
                     })

# compute stacking weights via loo
W <- loo_model_weights(log_lik_list,
                       method = 'stacking',
                       r_eff_list = r_eff_list,
                       optim_control = list(reltol = 1e-10))

### GENERATE PREDICTIONS ###

# generate test or holdout data set
Rtest = 50
Ttest = 10

test_vars <- list(R = Rtest, T = Ttest, A = nalts, L = nlvls, C = ncovs)
test_data <- stan("./STAN/generate_data.stan",
                  data = test_vars,
                  iter = 1,
                  chains = 1,
                  warmup = 0,
                  refresh = 1,
                  seed = 1750532,
                  algorithm = 'Fixed_param')

testY <- extract(test_data, permuted = TRUE)$Y
testX <- extract(test_data, permuted = TRUE)$X

# get model coefficients
B_list <- lapply(fit_list,
                 function(x) {
                   extract(x, pars='B')$B
                 })

prediction_fit_list <- list()
for (k in 1:K) {
  # input data for the prediction models
  testing_vars <- list(A = nalts,
                       L = nlvls,
                       I = niter,
                       R = nresp,
                       Rtest = Rtest,
                       Ttest = Ttest,
                       B = as.array(B_list[[k]]),
                       Xtest = as.array(testX[1,,,,]))

  fit <- stan("./STAN/base_model_predictions.stan",
              data = testing_vars,
              iter = 100,
              chains = 1,
              warmup = 0,
              algorithm = 'Fixed_param')

  prediction_fit_list[[k]] <- fit
}

weights <- as.vector(W)
Y_count_list <- lapply(prediction_fit_list,
                       function(x) {
                         Yc <- extract(x, pars='Y_count')$Y_count
                         colSums(Yc, dims=1)
                       })

for (k in 1:K) {
  Y_count_list[[k]] <- weights[[k]] * Y_count_list[[k]]
}

total_Y_count <- Reduce("+", Y_count_list)

print(dim(total_Y_count))

hit_count <- 0

for (r in 1:Rtest) {
  for (t in 1:Ttest) {
    Y_predicted <- which.max(total_Y_count[r, t, ])
    if (Y_predicted == testY[1, r, t]) {
      hit_count <- hit_count + 1
    }
  }
}

print(hit_count/(Rtest*Ttest))

### END ###
