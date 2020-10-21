# Preamble ----------------------------------------------------------------
# Load packages.
library(tidyverse)
library(rstan)
library(bayesplot)
library(tidybayes)
library(loo)

# Set Stan options.
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Load data.

ind_non <- 1 # Indicates no pathologies.
ind_scr <- 0 # Indicates screening.
ind_ana <- 0 # Indicates attribute non-attendance.
ind_sna <- 0 # Indicates screening and attribute non-attendance.

# Load and format data.
if (ind_non == 1) {
  # Matrix of observed choices.
  Y <- read_csv(here::here("Data", "01_PathologyNone", "Y.csv")) %>%
    select(-resp) %>%
    as.matrix()

  # Array of experimental designs per choice task.
  X_raw <- read_csv(here::here("Data", "01_PathologyNone", "X.csv"))
}

# Restructure X_raw into X.
X <- array(NA, dim = c(max(X_raw$resp), max(X_raw$task), max(X_raw$alt), ncol(X_raw) - 3))
for (n in 1:max(X_raw$resp)) {
  for (t in 1:max(X_raw$task)) {
    X[n,t,,] <- X_raw %>%
      filter(resp == n, task == t) %>%
      select(contains("l_")) %>%
      as.matrix()
  }
}

# Matrix of respondent-level covariates.
Z <- rep(1, nrow(Y)) %>%
  as.matrix()

# load(here::here("Data", "R05_Zerorez", "design.RData"))
# Y <- data$train_Y
# X <- data$train_X
# Z <- matrix(rep(1, nrow(Y)), ncol = 1)
# mat_ana <- data$mat_ana

mat_ana <- mat_ana[,1:dim(X)[4]]

# Ensemble Calibration ----------------------------------------------------
K <- nrow(mat_ana)
ensemble_fit <- vector(mode = "list", length = K)
for (k in 1:K) {
  data <- list(
    R = dim(X)[1],    # Number of respondents.
    S = dim(X)[2],    # Number of choice tasks.
    A = dim(X)[3],    # Number of choice alternatives.
    I = dim(X)[4],    # Number of observation-level covariates.
    J = ncol(Z),      # Number of population-level covariates.
    K = K,            # Number of members of the ensemble.
    k = k,            # Ensemble member number.
    
    Gamma_mean = 0,   # Mean of population-level means.
    Gamma_scale = 5,  # Scale of population-level means.
    Omega_shape = 2,  # Shape of population-level scale.
    tau_mean = 0,     # Mean of population-level scale.
    tau_scale = 5,    # Scale of population-level scale.
    
    Y = Y,            # Matrix of observations.
    X = X,            # Array of observation-level covariates.
    Z = Z,            # Matrix of population-level covariates.
    mat_ana = mat_ana # Clever randomization matrix.
  )
  
  fit = stan(
    here::here("Code", "Source", "hmnl_ensemble.stan"),
    data = data,
    seed = 42
  )
  
  # Calibrate the model using variational inference.
  ensemble_fit[[k]] <- vb(
    stan_model(here::here("Code", "Source", "hmnl_ensemble.stan")),
    data = data,
    seed = 42
  )
}

# Save ensemble output.
write_rds(
  ensemble_fit,
  here::here("Output", "ensemble_fit_01.rds")
)



# Check observation model trace plots.
beta_string <- str_c("Beta[", 1:data$R, ",", 1, "]")
for (i in 2:data$I) {
  beta_temp <- str_c("Beta[", 1:data$R, ",", i, "]")
  beta_string <- c(beta_string, beta_temp)
}
beta_string_filter <- beta_string[1:(data$I * 10)]

# Beta.
fit %>%
  # ...filter so we just have a subset of respondents to confirm fixed values
  # (see tidybayes)
  mcmc_trace(
    pars = beta_string_filter,
    n_warmup = 500,
    facet_args = list(
      ncol = 12,
      labeller = label_parsed
    )
  )

ggsave(
  "mcmc_trace-beta.png",
  path = here::here("Figures"),
  width = 20, height = 30, units = "in"
)

which(mat_ana[1,]==1)

