# Preamble ----------------------------------------------------------------
# Load libraries.
library(tidyverse)
library(mvtnorm)
library(rstan)
library(bayesplot)

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
  beta <- rmvnorm(1, mean = Gamma %*% z_resp, sigma = Vbeta)
    
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
Data <- list(nresp = nresp, nscns = nscns, nalts = nalts, nlvls = nlvls, ncovs = ncovs,
             Y = Y, X = X, Z = Z, Gamma = Gamma, Vbeta = Vbeta, Beta = Beta)

# Run model.
out_test <- stan(
  file = "../MODELS/HBMNL_1.1.stan", data = Data, iter = 500, chains = 4
)

# Trace plots.
plot(out_test, plotfun = "trace", pars = "Gamma")
# plot(out_test, plotfun = "trace", pars = "tau")
# plot(out_test, plotfun = "trace", pars = "Omega")

# Interval plots.
plot(out_test, pars = "Gamma")
# plot(out_test, pars = c("Gamma", "tau", "Omega"))
# plot(out_test, pars = c("Gamma", "Vbeta"))

# Numeric summary.
summary(out_test, pars=c("Gamma"))$summary

