# This code implements a hierarchical MNL choice model
# with a multivariate normal distribution of heterogeneity via Stan.

# Preamble ----------------------------------------------------------------
# Load libraries.
library(parallel)
library(loo)
library(rstan)

# Set Stan options.
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Parameter Recovery ------------------------------------------------------
nresp <- 100    # Number of respondents.
nscns <- 10     # Number of choice scenarios for each respondent.
nalts <- 4      # Number of alternatives in each choice scenario.
nlvls <- 12     # Number of attribute levels for each alternative.
ncovs <- 1      # Number of covariates for each respondent.

# True Gamma (ncovs x nlvls) values and Vbeta (nlvls x nlvls) values.
Gamma <- matrix(runif(ncovs * nlvls, -3, 4), nrow = nlvls, ncol = ncovs)
Vbeta <- diag(nlvls) + .5 * matrix(1, nrow = nlvls, ncol = 1) %*% t(matrix(1, nrow = nlvls, ncol = 1))


# Generate data.
Y <- matrix(nrow = nresp, ncol = nscns)
X <- array(dim = c(nresp, nscns, nalts, nlvls))
Z <- matrix(nrow = ncovs, ncol = nresp)
Beta <- matrix(nrow = nlvls, ncol = nresp)
for (resp in 1:nresp) {
  # Generate covariates for the distribution of heterogeneity.
  z_resp <- 1
  if (ncovs > 1) z_resp <- c(z_resp, round(runif(ncovs - 1)))
  
  # Generate individual-level betas.
  beta <- Gamma %*% z_resp + chol(Vbeta) %*% rnorm(nlvls)
    
  # Compute the latent utility a scenario at a time.
  for (scn in 1:nscns) {
    # Generate the design matrix for the given scenario.
    X_scn <- matrix(round(runif(nalts * nlvls)), nrow = nalts, ncol = nlvls)
    
    # Compute and the latent utility for each alternative and find the max.
    U_scn <- X_scn %*% as.vector(beta) + matrix((-log(-log(runif(nalts)))), ncol = 1)
    
    # Save each scenario's data.
    Y[resp, scn] <- which(U_scn == max(U_scn))
    X[resp, scn, , ] <- X_scn
  }
  
  # Save out each respondent's data.
  Z[, resp] <- as.vector(z_resp)
  Beta[, resp] <- as.vector(beta)
}

# Save data in a list.
Data <- list(J = nresp, S = nscns, C = nalts, K = nlvls, G = ncovs,
             Y = Y, X = X, Z = t(Z), Gamma = Gamma, Vbeta = Vbeta,
             Beta = Beta, w = rbinom(nlvls, 1, .5))

log_lik_list = log_lik_list_temp <- list()
K <- 3
for (k in 1:K){
    # Fit the k-th model with Stan
    fit <- stan("./MODELS/HBMNL_02.stan", data = Data, chains = 2, iter = 800, cores = 4)

    log_lik_list[[k]] <- extract(fit)[["log_lik"]]
    log_lik_list_temp[[k]] <- matrix(
      NA, 
      nrow = dim(log_lik_list[[k]])[1], 
      ncol = dim(log_lik_list[[k]])[2] * dim(log_lik_list[[k]])[3]
    )
    for (i in 1:dim(log_lik_list[[k]])[1]) {
      log_lik_list_temp[[k]][i,] <- as.vector(log_lik_list[[k]][1,,])
    }
}
log_lik_list <- log_lik_list_temp

model_weights <- loo_model_weights(log_lik_list, method="stacking")
print(model_weights)

