<<<<<<< HEAD
predictive_fit_ensemble = function(ensemble_weights, ensemble_fit, test_X, test_Y, Z){
=======
predictive_fit_ensemble = function(ensemble_weights, ensemble_fit, test_X, test_Y){
>>>>>>> e360f7c9931ecc39cc0d14ea388bdea47f5ad2e9
  # Compute the hit rate, hit prob, and loo metrics for the ensemble model.
  #   ensemble_weights - estimated weights for each of the models
  #   ensemble_fit - ensemble output with log_lik, betadraws, gammas, and Omegas for each model
  #   test_Y - choices (hold-out sample)
  #   test_X - design matrices (hold-out sample)
  #   Z - matrix of covariates

  nens <- length(ensemble_fit)
  ndraw <- length(ensemble_fit[[1]]$Beta[,1,1])         # Number of draws
  nresp <- length(test_Y[,1])           # Number of respondents
  nscns <- length(test_X[1, ,1,1])      # Number of choice tasks
  nalts <- length(test_X[1,1, ,1])      # Number of alternatives 
  nlvls <- length(test_X[1,1,1, ])      # Number of att levels
  if(is.null(Z)){Z <- matrix(double(nresp)+1, nc = 1)}
  ncov <- ncol(Z)  # Number of covariates
  
  #weight log_lik for each model to get log_lik for ensemble
<<<<<<< HEAD
  LLarray_ens = array(0, dim(ensemble_fit[[1]]$log_lik))
=======
  LLarray_ens = array(0,dim(ensemble_draws[[1]]$log_lik))
>>>>>>> e360f7c9931ecc39cc0d14ea388bdea47f5ad2e9
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
  probs_ens_postbetas = matrix(0, nr = nalts , nc = nresp*nscns*ndraw)
  probs_ens_Zgamma = matrix(0, nr = nalts , nc = nresp*nscns*ndraw)
  probs_ens_rand_het = matrix(0, nr = nalts , nc = nresp*nscns*ndraw)
  
  #loop over models
  for(model in 1:nens){
    #get betas for 2 different hit rate calculations:
      #1) using mean of posterior betas as predicted part-worths
      Umat_postbetas = matrix(0, nr = nresp*nscns*nalts, nc=ndraw)
      #2) using mean of the distribution of heterogeneity as predicted part worths
      Umat_Zgamma = matrix(0, nr = nresp*nscns*nalts, nc=ndraw)
    
    #loop over respondents
    for(resp in 1:nresp){
      #pull betas
      betadraw <- ensemble_fit[[model]]$Beta
      #transpose dimensions of betadraw array
      betadraw <- aperm(betadraw, c(2,1,3))
      #convert to matrix for faster mean over respondents
      betadraw <- matrix(betadraw, nr= dim(betadraw)[1])
      #take mean over respondents to get mean posterior betas for this model
      post_betabar_mat<- matrix(colMeans(betadraw), nr=ndraw)
      #identify which betas are set to 0 in this model
      index <- double(nlvls)
      index[colSums(post_betabar_mat)!=0] <- 1
      
      #use upper level to get mean of dist of heterogeneity
      gammadraw=ensemble_fit[[model]]$Gamma
      #transpose dimensions of gammadraw array
      gammadraw <- aperm(gammadraw, c(2,1,3))
      #multiply by Z to get mean of dist of het
      Zgamma_mat <- matrix(Z[resp,]%*%gammadraw, nr=ndraw)
      #set Zgamma_mat column = 0 if ensemble ignores the level
      Zgamma_mat[,index==0] <- 0
      
      #get utility for each alternative using the three different types
      Umat_postbetas[((resp-1)*nalts*nscns+1):((resp)*nalts*nscns),] <- 
                           exp(test_X_stacked[((resp-1)*nalts*nscns+1):((resp)*nalts*nscns),]%*%
                                            t(post_betabar_mat))
      Umat_Zgamma[((resp-1)*nalts*nscns+1):((resp)*nalts*nscns),] <- 
        exp(test_X_stacked[((resp-1)*nalts*nscns+1):((resp)*nalts*nscns),]%*%
        t(Zgamma_mat))
    }
      
    #find probabilities for each task, resp and draw for Postbetas
    Umat_postbetas <- matrix(Umat_postbetas, nr = nalts) 
    sums <- t(matrix(rep(colSums(Umat_postbetas),nalts), nc=nalts))
    #combine with other model probs weight by ensemble weights 
    probs_ens_postbetas <- probs_ens_postbetas + 
                              (Umat_postbetas/sums) *ensemble_weights[model]
      
    #find probabilities for each task, resp and draw for Zgamma
    Umat_Zgamma <- matrix(Umat_Zgamma, nr = nalts) 
    sums <- t(matrix(rep(colSums(Umat_Zgamma),nalts), nc=nalts))
    #combine with other model probs weight by ensemble weights 
    probs_ens_Zgamma <- probs_ens_Zgamma + 
      (Umat_Zgamma/sums) *ensemble_weights[model]
  }
  
  #find location of highest prob postbetas
  locs_postbetas <- apply(probs_ens_postbetas,2,which.max)
  
  #find location of highest prob Zgamma
  locs_Zgamma <- apply(probs_ens_Zgamma,2,which.max)
      
  #calculate hits postbetas
  hits_postbetas <- double(nresp*nscns*ndraw)
  hits_postbetas[locs_postbetas==rep(test_Y_stacked,ndraw)] <- 1
  
  #calculate hit probs postbetas
  hit_probs_postbetas <- colSums(probs_ens_postbetas*diag(nalts)[,rep(test_Y_stacked,ndraw)])
  
  #calculate hits Zgamma
  hits_Zgamma <- double(nresp*nscns*ndraw)
  hits_Zgamma[locs_Zgamma==rep(test_Y_stacked,ndraw)] <- 1
  
  #calculate hit probs Zgamma
  hit_probs_Zgamma <- colSums(probs_ens_Zgamma*diag(nalts)[,rep(test_Y_stacked,ndraw)])
  
  hitrates=c(mean(hits_postbetas),mean(hits_Zgamma))
  hitprobs=c(mean(hit_probs_postbetas), mean(hit_probs_Zgamma))

  return(list(hit_prob=hitprobs, hit_rate=hitrates, loo_fit=loo_fit_ens))
}



