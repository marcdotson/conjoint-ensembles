predictive_fit_ensemble = function(ensemble_weights, ensemble_fit, test_X, test_Y){
  # Compute the hit rate, hit prob, and loo metrics for the ensemble model.
  #   ensemble_weights - estimated weights for each of the models
  #   ensemble_fit - ensemble output with log_lik and betadraws for each model
  #   test_Y - choices (hold-out sample)
  #   test_X - design matrices (hold-out sample)
  
  n_ens <- length(ensemble_fit)
  ndraw <- length(ensemble_fit[[1]]$Beta[,1,1])         # Number of draws
  nresp <- length(test_Y[,1])           # Number of respondents
  nscns <- length(test_X[1, ,1,1])      # Number of choice tasks
  nalts <- length(test_X[1,1, ,1])      # Number of alternatives 
  nlvls <- length(test_X[1,1,1, ])      # Number of att levels
  
  #weight log_lik for each model to get log_lik for ensemble
  LLarray_ens = array(0,dims)
  for(k in 1:n_ens){
    #extract log_lik array from each stanfit object
    LLarray_ens <- LLarray_ens + weights[k]*ensemble_fit[[k]]$log_lik
  }  
  
  #get effective sample size
  r_eff_ens <- loo::relative_eff(x = exp(LLarray_ens), cores = cores)
  
  #apply PSIS via loo to ensemble likelihoods (loo fit metrics)
  loo_fit_ens <- loo::loo.array(LLarray, r_eff = r_eff,
                                cores = cores, save_psis = FALSE)
  
  #stack resps and scns to avoid loops (this needs changed if using hold out tasks)
  test_X_stacked <- NULL
  for(resp in 1:nresp){
    for(scn in 1:nscns){
      test_X_stacked <- rbind(test_X_stacked,test_X[resp,scn,,])
    }
  }
  
  #stack scn choices to avoid loops
  test_Y_stacked <- matrix(t(test_Y),nr=1)
  
  #loop over ensemble models to calculate individual hit rates
  hit_rate_vec=double(n_ens)
  hit_prob_vec=double(n_ens)
  for(n in 1:n_ens){
    betadraw <- ensemble_fit[[n]]$Beta
    
    #average betadraws over in-sample respondents for hold-out sample (no upper level)
    betadraw_avg <- matrix(double(ndraw*nlvls), nc=nlvls)
    for(draw in 1:ndraw){
      betadraw_avg[draw,] <- colMeans(betadraw[draw,,])
    }
    
    #find probabilities for each
    Umat <- exp(test_X_stacked%*%t(betadraw_avg))
    Umat <- matrix(Umat, nr = 3) 
    sums <- t(matrix(rep(colSums(Umat),3), nc=3))
    probs <- Umat/sums
    
    #find location of highest prob
    locs <- apply(probs,2,which.max)
    
    #identify hitrate and store
    hits <- double(length(locs))
    hits[ locs== rep(test_Y_stacked,ndraw) ] <- 1
    hit_rate_vec[n]=mean(hits)
    
    #identify hitprobs and store
    prob_select <- colSums(probs*diag(nalts)[,rep(test_Y_stacked,ndraw)])
    hit_prob_vec[n]=mean(prob_select)
  }
  
  hit_rate_ens <- hit_rate_vec%*%ensemble_weights
  hit_prob_ens <- hit_prob_vec%*%ensemble_weights
  return(list(hit_prob=hit_prob_ens, hit_rate=hit_rate_ens, loo_fit=loo_fit_ens))
}
