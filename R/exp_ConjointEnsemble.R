# R script to run a pathology treatment ensemble # Preamble ----------------------------------------------------------------
# Load libraries.
library(rstan)

# For fold in fold_list:
#   combine other 4 folds as train set
#   for each base model:
#       fit base model to training fold
#       make predictions on test fold and save in training fold
# Fit each base model to full training data set
# make predictions on test data set and save in test meta
# Fit new model on training data with base model predictions included as features
# use new model to make predictions on test meta

# Set Stan options.
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

nresp = 100
ntask_test = 20
ntask = 10
nalts = 4
nlvls = 12
ncovs = 1
niter = 300
nchains = 2
treedepth = 3
random_seed = 1750532
model_name <- "./STAN/mnl_vanilla.stan"
Ttest <- 10
K <- 2

# input data for generate_data.stan
train_vars <- list(R = nresp, T = ntask_test, A = nalts, L = nlvls, C = ncovs)

# Generate training data set
train_data <- stan("./STAN/generate_data.stan",
                   data = train_vars,
                   iter = 1,
                   chains = 1,
                   warmup = 0,
                   refresh = 1,
                   seed = random_seed,
                   algorithm = 'Fixed_param')
                  
# Extract training data from fit
simX <- extract(train_data, permuted = TRUE)$X
simY <- extract(train_data, permuted = TRUE)$Y
simZ <- extract(train_data, permuted = TRUE)$Z

testX <- as.array(simX[1,,11:20,,])
testY <- as.array(simY[1,,11:20])

trainX <- as.array(simX[1,,1:10,,])
trainY <- as.array(simY[1,,1:10])

# Fit base models to training data
fit_list <- list()
for (k in 1:K) {
  t_ <- 1+(k-1)*5
  T_ <- 5+(k-1)*5
  # input data for the Stan base models
  training_vars <- list(R = nresp,
                        T = 5,
                        A = nalts,
                        L = nlvls,
                        C = ncovs,
                        X = trainX[,t_:T_,,],
                        Y = trainY[,t_:T_],
                        Z = as.matrix(simZ[1,,]))

  # Fit the k-th model with Stan
  fit <- stan(model_name,
              data = training_vars,
              iter=niter,
              chains = nchains,
              control = list(adapt_delta = .9, max_treedepth = treedepth))
  fit_list[[k]] <- fit
}

### GENERATE PREDICTIONS ###

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
                       I = (nchains/2)*niter,
                       R = nresp,
                       Ttest = Ttest,
                       B = as.array(B_list[[k]]),
                       Xtest = testX)

  fit <- stan("./STAN/model_predictions.stan",
              data = testing_vars,
              iter = 100,
              chains = 1,
              warmup = 0,
              algorithm = 'Fixed_param')

  prediction_fit_list[[k]] <- fit
}

Y_pp_list <- lapply(prediction_fit_list,
                    function(x) {
                      Yc <- extract(x, pars='Y_pp')$Y_count
                      colSums(Yc, dims=1)
                    })

total_Y_count <- Reduce("+", Y_count_list)


print(dim(total_Y_count))

hit_count <- 0

for (r in 1:nresp) {
  for (t in 1:Ttest) {
    Y_predicted <- which.max(total_Y_count[r, t, ])
    if (Y_predicted == testY[r, t]) {
      hit_count <- hit_count + 1
    }
  }
}

print(hit_count/(nresp*Ttest))

### END ###
