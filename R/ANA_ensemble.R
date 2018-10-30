# This code implements a hierarchical MNL choice model
# with a multivariate normal distribution of heterogeneity via Stan.

# Preamble ----------------------------------------------------------------
# Load libraries.
Sys.setenv(USE_CXX14 = 1)
library(parallel)
library(loo)
library(rstan)

# Set Stan options.
#rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

pathology <- function(beta, nlvls) {
  beta <- beta*rbinom(nlvls, 1, .5);
  return(beta)
}

# Parameter Recovery ------------------------------------------------------
nresp <- 110    # Number of respondents.
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
Z <- matrix(nrow = nresp, ncol = ncovs)
Beta <- matrix(nrow = nlvls, ncol = nresp)

for (resp in 1:nresp) {
  # Generate covariates for the distribution of heterogeneity.
  z_resp <- 1
  if (ncovs > 1) z_resp <- c(z_resp, round(runif(ncovs - 1)))
  
  # Generate individual-level betas.
  beta <- Gamma %*% z_resp + chol(Vbeta) %*% rnorm(nlvls)
  # implement Pathology
  beta <- pathology(beta, nlvls)
    
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
  Z[resp, ] <- as.vector(z_resp)
  Beta[, resp] <- as.vector(beta)
}

# Save data in a list.
Data <- list(J = 100, S = nscns, C = nalts, K = nlvls, G = ncovs,
             Y = Y[1:100, ], X = X[1:100, , , ], Z = as.matrix(Z[1:100, 1]), Gamma = Gamma, Vbeta = Vbeta,
             Beta = Beta[, 1:100], w = rbinom(nlvls, 1, .5),
             Ytest = Y[100:nresp, ], Xtest = X[100:nresp, , , ], Betatest = Beta[, 100:nresp])

log_lik_list <- list()
K <- 3
for (k in 1:K){
    # Fit the k-th model with Stan
    fit <- stan("./MODELS/HBMNL_ana2.stan", data = Data, chains = 2, iter = 300, cores = 4)
    log_lik_list[[k]] <- extract_log_lik(fit)
}
model_weights <- loo_model_weights(log_lik_list, method="stacking")
print(model_weights)
