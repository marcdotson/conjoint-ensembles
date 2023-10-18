#********************************************************
#
#  ensemble_weights
#  Function for the calculation of ensemble model weights
#
#********************************************************

#Input
# X - List of Stanfit objects for each model 
#     in the ensemble.  Must include log_lik

#Output
# weights - vector of model weights

#function to calculate weights
ensemble_weights <- function(X, cores){
  PSIS_list <- NULL  #list of PSIS loo objects
  K <- length(X)
  for(k in 1:K){
    #extract log_lik array from each stanfit object
    LLarray <- loo::extract_log_lik(stanfit = X[[k]],
                                  parameter_name = 'log_lik',
                                  merge_chains = FALSE)
    #get relative effective sample size for array
    r_eff <- loo::relative_eff(x = exp(LLarray), cores = cores)
    #apply PSIS via loo to array and save
    PSIS_list[[k]] = loo::loo.array(LLarray, r_eff = r_eff,
                                    cores = cores, save_psis = FALSE)
  }  
  set.seed(42)
  #calculate weights
  weights = loo::loo_model_weights(x = PSIS_list, method = "stacking", 
                              optim_method = "BFGS", optim_control = list(reltol=1e-10),
                              r_eff_list = r_eff_list, cores = cores, )
  return(weights)
}

