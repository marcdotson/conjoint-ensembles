predictive_fit_ensemble = function(ensemble_weights, ensemble_fit, test_X, test_Y){
  # Compute the hit rate, hit prob, and loo metrics for the ensemble model.
  #   ensemble_weights - estimated weights for each of the models
  #   ensemble_fit - ensemble output with log_lik and betadraws for each model
  #   test_Y - choices (hold-out sample)
  #   test_X - design matrices (hold-out sample)
  
  nens <- length(ensemble_fit)
  ndraw <- length(ensemble_fit[[1]]$Beta[,1,1])         # Number of draws
  nresp <- length(test_Y[,1])           # Number of respondents
  nscns <- length(test_X[1, ,1,1])      # Number of choice tasks
  nalts <- length(test_X[1,1, ,1])      # Number of alternatives 
  nlvls <- length(test_X[1,1,1, ])      # Number of att levels
  
  #weight log_lik for each model to get log_lik for ensemble
  LLarray_ens = array(0,dim(ensemble_draws[[1]]$log_lik))
  for(k in 1:nens){
    #extract log_lik array from each stanfit object
    LLarray_ens <- LLarray_ens + ensemble_weights[k]*ensemble_fit[[k]]$log_lik
  }  
  
  #get effective sample size
  cores <- parallel::detectCores()
  r_eff_ens <- loo::relative_eff(x = exp(LLarray_ens), cores = cores)
  
  #apply PSIS via loo to ensemble likelihoods (loo fit metrics)
  loo_fit_ens <- loo::loo.array(LLarray_ens, r_eff = r_eff_ens,
                                cores = cores, save_psis = FALSE)
  
  #stack resps and scns to avoid loops (this needs changed if using hold out tasks)
  test_X_stacked <- NULL
  for(resp in 1:nresp){
    for(scn in 1:nscns){
      test_X_stacked <- rbind(test_X_stacked,test_X[resp,scn,,])
    }
  }
  
  #stack scn choices to avoid loops
  test_Y_stacked <- matrix(t(test_Y),nc=1)
  
  #loop over ensemble models to calculate individual hit rates
  hit_rate_vec=double(ndraw)
  hit_prob_vec=double(ndraw)
  #loop over draws
  for(draw in 1:ndraw){
    #pull the posterior mean betas from each model into a matrix for this draw
    betabar_mat <- matrix(double(nlvls*nens), nc=nlvls)
    for(n in 1:nens){
      betadraw <- ensemble_fit[[n]]$Beta
      #average betadraws over in-sample respondents for hold-out sample (no upper level)
      betabar_mat[n,] <- colMeans(betadraw[draw,,])
    }
    
    #find probabilities for each
    Umat <- exp(test_X_stacked%*%t(betabar_mat))
    Umat <- matrix(Umat, nr = nalts) 
    sums <- t(matrix(rep(colSums(Umat),nalts), nc=nalts))
    probs <- Umat/sums
    
    #stack probs to weight by ensemble weights
    probs <- matrix(probs, nc=nens)
    
    #get probs for ensemble and reorganize by scn
    probs_ens <- probs%*%ensemble_weights
    probs_ens <- matrix(probs_ens, nr=nalts)
    
    #find location of highest prob
    locs <- apply(probs_ens,2,which.max)
    
    #identify hitrate and store
    hits <- double(length(locs))
    hits[ locs == test_Y_stacked ] <- 1
    hit_rate_vec[draw]=mean(hits)
    
    #identify hitprobs and store
    prob_select <- colSums(probs_ens*diag(nalts)[,test_Y_stacked])
    hit_prob_vec[draw] <- mean(prob_select)
  }
  
  hit_rate_ens <- mean(hit_rate_vec)
  hit_prob_ens <- mean(hit_prob_vec)
  return(list(hit_prob=hit_prob_ens, hit_rate=hit_rate_ens, loo_fit=loo_fit_ens))
}
