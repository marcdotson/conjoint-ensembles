predictive_fit_hmnl <- function(hmnl_draws, test_X, test_Y, test_Z) {
  # Compute the hit rate for the indicated model.
  #   hmnl_draws - hmnl output with Beta, Gamma, Omega, and tau draws
  #   test_Y - choices (hold-out sample)
  #   test_X - design matrices (hold-out sample)
  #   test_Z - matrix of covariates

  nresp <- dim(test_X)[1] # Number of respondents.
  nscns <- dim(test_X)[2] # Number of choice tasks.
  
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
  
  # Get the mean of distribution of heterogeneity.
  # Gamma_mean <- apply(hmnl_draws$Gamma, c(2,3), mean)
  Gamma_mean <- hmnl_draws |> 
    subset_draws(variable = "Gamma") |> 
    summarize_draws("mean") |> 
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
  
  # Return the average hit rate and hit probability.
  return(list(hit_prob = mean(probs), hit_rate = mean(hits)))
}
