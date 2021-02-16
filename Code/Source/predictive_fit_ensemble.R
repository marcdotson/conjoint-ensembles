predictive_fit_ensemble = function(indices, ensemble_weights, ensemble_draws, 
                                   test_X, test_Y, mat_ana, mat_screen, Z){
  # Compute the hit rate, hit prob, and loo metrics for the ensemble model.
  #   ensemble_weights - estimated weights for each of the models
  #   ensemble_fit - ensemble output with log_lik, betadraws, gammas, and Omegas for each model
  #   test_Y - choices (hold-out sample)
  #   test_X - design matrices (hold-out sample)
  #   Z - matrix of covariates
  
  nens <- length(ensemble_draws)
  ndraw <- length(ensemble_draws[[1]]$log_lik[,1,1]) # Number of draws
  nens <- length(ensemble_draws)
  nresp <- length(test_Y[,1])           # Number of respondents
  nscns <- length(test_X[1, ,1,1])      # Number of choice tasks
  nalts <- length(test_X[1,1, ,1])      # Number of alternatives 
  nlvls <- length(test_X[1,1,1, ])      # Number of att levels
  if( is.null(Z) ) Z <- matrix(1, nr = nresp, nc = 1)
  
  #weight log_lik for each model to get log_lik for ensemble
  LLmat_ens = matrix(0, nr=ndraw , 
                     nc=exp(sum(log(dim(ensemble_draws[[k]]$log_lik))))/ndraw)
  loglik=ensemble_draws[[k]]$log_lik
  ndraw=dim(loglik)[1]
  for(k in 1:nens){
    #extract log_lik array from each stanfit object
    loglik=ensemble_draws[[k]]$log_lik
    LLmat <- matrix(loglik, nr=ndraw)
    LLmat_ens <- LLmat_ens + ensemble_weights[k]*LLmat
  }  
  
  #get effective sample size
  cores <- parallel::detectCores()
  r_eff_ens <- loo::relative_eff(x = exp(LLmat_ens), chain_id = double(ndraw)+1, cores = cores)
  
  #apply PSIS via loo to ensemble likelihoods (loo fit metrics)
  loo_fit_ens <- loo::loo.matrix(LLmat_ens, r_eff = r_eff_ens,
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
  
  #loop over ensemble models to calculate individual hit probs 
  probs_ens <- matrix(0, nr = nalts , nc = nresp*nscns)
  
  #loop over models
  for(model in 1:nens){
    #get betas for different hit rate calculations
    #using mean (over draws) of the mean of the post dist
    Umat <- matrix(0, nr = nresp*nscns*nalts)
    
    #get gammas
    meangammas=ensemble_draws[[model]]$Gamma
    
    #adjust due to pathology approach in ensembles
    index_ana <- mat_ana[model,]
    
    #loop over respondents
    for(resp in 1:nresp){
      #multiply by Z to get mean of dist of het for resp
      betas <- matrix(Z[resp,]%*%meangammas, nc=1)
      
      #set gammadraws_mat column = 0 if ensemble ignores the level
      betas[index_ana==1,] <- 0
      
      #get utility for each alternative
      Umat[((resp-1)*nalts*nscns+1):((resp)*nalts*nscns),] <- 
        exp(test_X_stacked[((resp-1)*nalts*nscns+1):((resp)*nalts*nscns),]%*%
              matrix(betas))
    }
    
    #find probabilities for each task and resp
    Umat_byscn <- matrix(Umat, nr = nalts) 
    sums <- t(matrix(rep(colSums(Umat_byscn),nalts), nc=nalts))
    #combine with other model probs weight by ensemble weights 
    probs_ens <- probs_ens + 
      (Umat_byscn/sums)*ensemble_weights[model]
  }

  #find location of highest prob meangammas
  locs <- apply(probs_ens,2,which.max)
  
  #calculate hits meangammas
  hits <- double(nresp*nscns)
  hits[locs==test_Y_stacked] <- 1
  
  #calculate hit probs
  hit_probs <- colSums(probs_ens*diag(nalts)[,test_Y_stacked])
  
  return(list(hit_rate = mean(hits), hit_prob = mean(hit_probs), loo_fit=loo_fit_ens))
}

