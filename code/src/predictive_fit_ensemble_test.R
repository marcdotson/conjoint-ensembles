predictive_fit_ensemble <- function(indices, ensemble_weights, ensemble_draws, test_X, test_Y, mat_ana, mat_screen, test_Z, ensemble_fit) {
  # Compute the hit rate, hit prob, and loo metrics for the ensemble model.
  #   ensemble_weights - estimated weights for each of the models
  #   ensemble_draws - ensemble output with log_lik, betadraws, gammas, and Omegas for each model
  #   test_Y - choices (hold-out sample)
  #   test_X - design matrices (hold-out sample)
  #   test_Z - matrix of covariates
  
  nmemb <- length(ensemble_draws) # Number of ensemble members.
  nresp <- dim(test_X)[1]         # Number of respondents.
  nscns <- dim(test_X)[2]         # Number of choice tasks.
  
  # Log-likelihood function for the MNL.
  ll_mnl = function (beta, y, X) {
    nvars = ncol(X)        # Number of attribute levels.
    nscns = length(y)      # Number of choice tasks.
    nalts = nrow(X)/nscns  # Number of alternatives.
    
    # Compute Xbeta across all choice tasks.
    Xbeta = matrix(exp(X%*%beta),byrow=TRUE,ncol=nalts)
    
    # Numerator: Xbeta values associated with each choice.
    choices = cbind(c(1:nscns),y)
    numerator = Xbeta[choices]
    
    # Denominator: Xbeta values associated with each task.
    iota = c(rep(1,nalts))
    denominator = Xbeta%*%iota
    
    # Return the logit summed across choice tasks.
    return(sum(log(numerator) - log(denominator)))
  }
  
  # Compute LOO and predict Y.
  memb_loo <- NULL
  memb_hits <- NULL
  memb_probs <- NULL
  for (memb in 1:nmemb) {
    # Get the mean of distribution of heterogeneity.
    # Gamma_mean <- ensemble_draws[[memb]]$Gamma
    Gamma_mean <- ensemble_fit$ensemble_draws[[1]]$Gamma |> 
      select("mean") |> 
      as.matrix()
    
    # Predict Y.
    hits <- NULL
    probs <- NULL
    for (resp in 1:nresp) {
      for (scn in 1:nscns) {
        # Pull the relevant Y and X.
        Y_scn <- test_Y[resp, scn]
        X_scn <- test_X[resp, scn,,]
        
        # Compute hit and probability.
        # hits <- c(hits, Y_scn == which.max(X_scn %*% t(Gamma_mean)))
        # probs <- c(probs, exp(ll_mnl(t(Gamma_mean), Y_scn, X_scn)))
        hits <- c(hits, Y_scn == which.max(X_scn %*% Gamma_mean))
        probs <- c(probs, exp(ll_mnl(Gamma_mean, Y_scn, X_scn)))
      }
    }
    
    # Compute LOO, the average hit rate, and hit probability.
    # memb_loo <- c(memb_loo, loo(ensemble_fit[[memb]])$elpd_loo)
    memb_loo <- c(memb_loo, 0)
    memb_hits <- c(memb_hits, mean(hits))
    memb_probs <- c(memb_probs, mean(probs))
  }
  
  # Return the weighted average hit rate and hit probability.
  return(list(
    hit_rate = ensemble_weights %*% memb_hits, 
    hit_prob = ensemble_weights %*% memb_probs, 
    loo_fit = ensemble_weights %*% memb_loo
  ))
}
