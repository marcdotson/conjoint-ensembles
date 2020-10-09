# START WITH EXISTING DATA; HOW WAS IT SIMULATED?
# NO NEED TO SPLIT INTO TRAIN/TEST WHEN USING LOO.
# IMPLEMENT MODIFICATION TO SIMULATE FOR PATHOLOGIES.

# Preamble ----------------------------------------------------------------
# Load libraries.
library(tidyverse)
library(mvtnorm)
library(rstan)
library(bayesplot)

# Set Stan options.
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Set the seed.
set.seed(42)

ind_non <- 0 # Indicates no pathologies.
ind_scr <- 0 # Indicates screening.
ind_ana <- 0 # Indicates attribute non-attendance.
ind_sna <- 0 # Indicates screening and attribute non-attendance.

# Simulate Data -----------------------------------------------------------
N <- 500          # Number of respondents.
S <- 10           # Number of choice tasks per respondent.
P <- 4            # Number of alternatives per choice task.
L <- 12           # Number of (estimable) attribute levels.
C <- 1            # Number of covariates.

# True Theta (C x L) values and Sigma (L x L) values.
Theta <- matrix(runif(C * L, -3, 4), byrow = FALSE, nrow = C, ncol = L)
Sigma <- diag(L) + .5 * matrix(1, nrow = L, ncol = 1) %*% t(matrix(1, nrow = L, ncol = 1))

# Generate Y, X, Z, and Beta.
Y <- matrix(NA, nrow = N, ncol = S)
X <- array(NA, dim = c(N, S, P, L))
Z <- matrix(NA, nrow = N, ncol = C)
Beta <- matrix(NA, nrow = N, ncol = L)
Y_list = X_list <- NULL
for (n in 1:N) {
  # Generate covariates for the distribution of heterogeneity.
  Z[n, 1] <- 1
  if (C > 1) Z[n, -1] <- c(z_temp, round(runif(C - 1)))
  
  # Generate individual-level betas.
  Beta[n,] <- mvtnorm::rmvnorm(1, mean=t(Theta) %*% Z[n,], sigma = Sigma)
  
  # Compute the latent utility a scenario at a time.
  X_temp <- NULL
  for (s in 1:S) {
    # Generate or draw the design matrix for the given scenario.
    X[n, s,,] <- matrix(round(runif(P * L)), nrow = P, ncol = L)
    
    # Compute and the latent utility for each alternative and find the max.
    U <- X[n, s,,] %*% Beta[n,] + matrix((-log(-log(runif(P)))), ncol = 1)
    Y[n, s] <- which(U == max(U))
    
    # Save out each design matrix as a list.
    X_temp <- rbind(X_temp, X[n, s,,])
  }
  
  # Also save out each respondent's data as a list.
  Y_list[[n]] <- Y[n,]
  X_list[[n]] <- X_temp
}

sim_data <- list(
  Y = Y,
  X = X,
  Z = Z,
  Y_list = Y_list,
  X_list = X_list,
  Theta = Theta,
  Sigma = Sigma,
  Beta = Beta
)

readr::write_rds(sim_data, here::here("Data", "hmnl_sim_data.rds"))
