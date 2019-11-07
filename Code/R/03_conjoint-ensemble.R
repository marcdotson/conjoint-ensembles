# Preamble ----------------------------------------------------------------
# Load packages.
library(tidyverse)
library(rstan)
library(bayesplot)
library(tidybayes)

# Stan options.
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Load and format data.
Y <- read_csv(here::here("Data", "01_PathologyNone", "Y.csv")) %>% 
  select(-resp) %>% 
  as.matrix()

X_raw <- read_csv(here::here("Data", "01_PathologyNone", "X.csv"))
X <- array(NA, dim = c(max(X_raw$resp), max(X_raw$task), max(X_raw$alt), ncol(X_raw) - 3))
for (n in 1:max(X_raw$resp)) {
  for (t in 1:max(X_raw$task)) {
    X[n,t,,] <- X_raw %>% 
      filter(resp == n, task == t) %>% 
      select(contains("l_")) %>% 
      as.matrix()
  }
}

Z <- rep(1, nrow(Y)) %>% 
  as.matrix()

# Estimation --------------------------------------------------------------
# Specify the data for calibration in a list.
data <- list(
  N = dim(X)[1],           # number of respondents
  S = dim(X)[2],           # number of questions (unique inquiries)
  P = dim(X)[3],           # number of alternatives (choices) per question
  L = dim(X)[4],           # number of feature variables
  C = dim(Z)[2],           # number of respondent covariates (demographics, etc)
  
  Y = Y,                   # observed responses
  X = X,                   # matrix of attributes for each obs
  Z = Z,                   # vector of covariates for each respondent
  
  mu_loc = 0,              # location of the means of B
  mu_scale = 1,            # scale of the means of B
  alpha_loc = 0,           # location of the variance of B
  alpha_scale = 10,        # scale of the variance of B
  lkj_corr_shape = 5       # for correlation matrix hyperprior
)

# Calibrate the model.
fit <- stan(
  file = here::here("Code", "Stan","hmnl.stan"), # NO LONGER EXISTS WITH THIS PARAMETERIZATION...
  data = data,
  seed = 42
)

write_rds(fit, here::here("Output", "fit.RDS"))

# Model Checking ----------------------------------------------------------
# Check divergences.
library(bayesplot)
source(here::here("Code", "R", "stan_utility.R"))

check_div(fit)

as.matrix(fit) %>% 
  mcmc_scatter(
    pars = c("B[1,1]", "B[1,2]"), 
    np = nuts_params(fit),
    np_style = scatter_style_np(div_color = "green", div_alpha = 0.50)
  )

# Check the effective sample size.
check_n_eff(fit)

# Check the Rhat statistic.
check_rhat(fit)

# Check trace plots.
fit %>% 
  extract(
    inc_warmup = TRUE, 
    permuted = FALSE
  ) %>% 
  mcmc_trace(
    pars = c(
      "B[1,1]", "B[1,2]", "B[1,3]", "B[1,4]", "B[1,5]", "B[1,6]", 
      "B[1,7]", "B[1,8]", "B[1,9]", "B[1,10]", "B[1,11]", "B[1,12]"
    ),
    n_warmup = 1000,
    facet_args = list(nrow = 2, labeller = label_parsed)
  )

# Recover parameter values.
as.array(fit) %>% 
  mcmc_areas(pars = c("B[1,1]"))

# # Check divergences.
# library(bayesplot)
# source(here::here("Code", "R", "stan_utility.R"))
# 
# check_div(fit)
# 
# as.matrix(fit) %>% 
#   mcmc_scatter(
#     pars = c("B[1,1]", "B[1,2]"), 
#     np = nuts_params(fit),
#     np_style = scatter_style_np(div_color = "green", div_alpha = 0.50)
#   )
# 
# # Check the effective sample size.
# check_n_eff(fit)
# 
# # Check the Rhat statistic.
# check_rhat(fit)
# 
# # Check trace plots.
# fit %>% 
#   extract(
#     inc_warmup = TRUE, 
#     permuted = FALSE
#   ) %>% 
#   mcmc_trace(
#     pars = c(
#       "B[1,1]", "B[1,2]", "B[1,3]", "B[1,4]", "B[1,5]", "B[1,6]", 
#       "B[1,7]", "B[1,8]", "B[1,9]", "B[1,10]", "B[1,11]", "B[1,12]"
#     ),
#     n_warmup = 1000,
#     facet_args = list(nrow = 2, labeller = label_parsed)
#   )
# 
# # Recover parameter values.
# as.array(fit) %>% 
#   mcmc_areas(pars = c("B[1,1]"))

