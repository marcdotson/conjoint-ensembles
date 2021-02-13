predictive_fit_ensemble = function(indices, ensemble_weights, ensemble_draws, 
                                   test_X, test_Y, mat_ana, mat_screen, Z){
  # Compute the hit rate, hit prob, and loo metrics for the ensemble model.
  #   ensemble_weights - estimated weights for each of the models
  #   ensemble_fit - ensemble output with log_lik, betadraws, gammas, and Omegas for each model
  #   test_Y - choices (hold-out sample)
  #   test_X - design matrices (hold-out sample)
  #   Z - matrix of covariates
  
  nens <- length(ensemble_draws)
  ndraw <- length(ensemble_draws[[1]]$Gamma[,1,1]) # Number of draws
  nresp <- length(test_Y[,1])           # Number of respondents
  nscns <- length(test_X[1, ,1,1])      # Number of choice tasks
  nalts <- length(test_X[1,1, ,1])      # Number of alternatives 
  nlvls <- length(test_X[1,1,1, ])      # Number of att levels
  if( is.null(Z) ) Z <- matrix(1, nr = nresp, nc = 1)
  ncov <- ncol(Z)  # Number of covariates
  
  #weight log_lik for each model to get log_lik for ensemble
  LLarray_ens = array(0, dim(ensemble_draws[[1]]$log_lik))
  for(k in 1:nens){
    #extract log_lik array from each stanfit object
    LLarray_ens <- LLarray_ens + ensemble_weights[k]*ensemble_draws[[k]]$log_lik
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
  
  #loop over ensemble models to calculate individual hit probs two ways
  probs_ens_gammadraws <- matrix(0, nr = nalts , nc = nresp*nscns*ndraw)
  probs_ens_meangammas <- matrix(0, nr = nalts , nc = nresp*nscns)
  
  #loop over models
  for(model in 1:nens){
    #get betas for different hit rate calculations:
    #1) using mean of posterior dist of heterogeneity
    Umat_gammadraws = matrix(0, nr = nresp*nscns*nalts, nc=ndraw)
    #2) using mean (over draws) of the mean of the post dist
    Umat_meangammas = matrix(0, nr = nresp*nscns*nalts)
    
    #get gammas
    gammadraw=ensemble_draws[[model]]$Gamma
    
    #transpose dimensions of gammadraw array
    gammadraw <- aperm(gammadraw, c(2,1,3))
    
    #adjust due to pathology approach in ensembles
    index_ana <- mat_ana[model,]
    
    
    #loop over respondents
    for(resp in 1:nresp){
      #multiply by Z to get mean of dist of het for resp
      gammadraws_mat <- matrix(Z[resp,]%*%gammadraw, nr=ndraw)
      #set gammadraws_mat column = 0 if ensemble ignores the level
      gammadraws_mat[,index_ana==1] <- 0
      #get utility for each alternative
      Umat_gammadraws[((resp-1)*nalts*nscns+1):((resp)*nalts*nscns),] <- 
        exp(test_X_stacked[((resp-1)*nalts*nscns+1):((resp)*nalts*nscns),]%*%
              t(gammadraws_mat))
      Umat_meangammas[((resp-1)*nalts*nscns+1):((resp)*nalts*nscns),] <- 
        exp(test_X_stacked[((resp-1)*nalts*nscns+1):((resp)*nalts*nscns),]%*%
              matrix(colMeans(gammadraws_mat)))
    }
    
    #find probabilities for each task, resp and draw for gammadraws
    Umat_gammadraws_byscn <- matrix(Umat_gammadraws, nr = nalts) 
    sums <- t(matrix(rep(colSums(Umat_gammadraws_byscn),nalts), nc=nalts))
    #combine with other model probs weight by ensemble weights 
    probs_ens_gammadraws <- probs_ens_gammadraws + 
      (Umat_gammadraws_byscn/sums)*ensemble_weights[model]
    
    #find probabilities for each task, resp and draw for gammadraws
    Umat_meangammas_byscn <- matrix(Umat_meangammas, nr = nalts) 
    sums <- t(matrix(rep(colSums(Umat_meangammas_byscn),nalts), nc=nalts))
    #combine with other model probs weight by ensemble weights 
    probs_ens_meangammas <- probs_ens_meangammas + 
      (Umat_meangammas_byscn/sums)*ensemble_weights[model]
  }
  
  #find location of highest prob gammadraws
  locs_gammadraws <- apply(probs_ens_gammadraws,2,which.max)
  #find location of highest prob meangammas
  locs_meangammas <- apply(probs_ens_meangammas,2,which.max)
  
  #calculate hits gammadraws
  hits_gammadraws <- double(nresp*nscns*ndraw)
  hits_gammadraws[locs_gammadraws==rep(test_Y_stacked,ndraw)] <- 1
  #calculate hits meangammas
  hits_meangammas <- double(nresp*nscns)
  hits_meangammas[locs_meangammas==test_Y_stacked] <- 1
  
  #calculate hit probs gammadraws
  hit_probs_gammadraws <- colSums(probs_ens_gammadraws*diag(nalts)[,rep(test_Y_stacked,ndraw)])
  #calculate hit probs meangammas
  hit_probs_meangammas <- colSums(probs_ens_meangammas*diag(nalts)[,test_Y_stacked])
  
  return(list(hit_prob_gammadraws=mean(hit_probs_gammadraws), hit_rate_gammadraws=mean(hits_gammadraws),
              hit_prob_meangammas=mean(hit_probs_meangammas), hit_rate_meangammas=mean(hits_meangammas),
              loo_fit=loo_fit_ens))
}

