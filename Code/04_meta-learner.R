# Preamble ----------------------------------------------------------------
# Load packages.
library(tidyverse)
library(rstan)
library(bayesplot)
library(tidybayes)
library(loo)

#read in data
ana_out <- readRDS("~/Desktop/reduced_fit-vb_ana.rds")

#extract list for ensembles
ensemble_fit <- ana_out$ensemble_fit

#variables
n_ens <- length(ensemble_fit) #number of ensembles
cores <- 4 # number of cores
dims <- dim(ensemble_fit[[1]]$log_lik) #dimensions of log_lik values


#create array of likelihoods with effective sample sizes
LooPSIS_list <- NULL  #list of PSIS loo objects
for(k in 1:n_ens){
  #extract log_lik array from each stanfit object
  LLarray <- ensemble_fit[[k]]$log_lik
  #get relative effective sample size for array
  r_eff <- loo::relative_eff(x = exp(LLarray), cores = cores)
  #apply PSIS via loo to array and save
  LooPSIS_list[[k]] = loo::loo.array(LLarray, r_eff = r_eff,
                                  cores = cores, save_psis = FALSE)
}  

#calculate weights
set.seed(22)
weights <- loo::loo_model_weights(x = LooPSIS_list, method = "stacking", 
                                 optim_method = "BFGS", optim_control = list(reltol=1e-10),
                                 r_eff_list = r_eff_list, cores = cores)
#prepare test data
test_Y <- ana_out$test_Y
test_X <- ana_out$test_X

#get functions for predictive fit
source(here::here("Code", "06_model-comparison.R"))
ana_fit <- predictive_fit_ensemble(ensemble_weights=weights, ensemble_fit=ensemble_fit, test_X, test_Y)
