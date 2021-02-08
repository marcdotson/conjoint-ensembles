predictive_fit_hmnl = function(hmnl_fit, test_X, test_Y, Z){
  # Compute the hit rate for the indicated model.
  #   hmnl_fit - hmnl output with log_lik, betadraws, gammadraws, and Omegadraws
  #   test_Y - choices (hold-out sample)
  #   test_X - design matrices (hold-out sample)
  #   Z - matrix of covariates

  ndraw <- length(hmnl_fit$Beta[,1,1]) # Number of draws
  nresp <- length(test_Y[,1])           # Number of respondents
  nscns <- length(test_X[1, ,1,1])      # Number of choice tasks
  nalts <- length(test_X[1,1, ,1])      # Number of alternatives 
  nlvls <- length(test_X[1,1,1, ])      # Number of att levels
  if(is.null(Z)){Z <- matrix(double(nresp)+1, nc = 1)}
  
  #stack resps and scns to avoid loops (this needs changed if using hold out tasks)
  test_X_stacked <- NULL
  for(resp in 1:nresp){
    for(scn in 1:nscns){
      test_X_stacked <- rbind(test_X_stacked,test_X[resp,scn,,])
    }
  }
  
  #stack scn choices to avoid loops
  test_Y_stacked <- matrix(t(test_Y),nc=1)
  
  #get utilities for 2 different hit rate calculations:
    #1) using mean of posterior betas as predicted part-worths
  Umat_postbetas = matrix(0, nr = nresp*nscns*nalts, nc=ndraw)
    #2) using mean of the distribution of heterogeneity as predicted part worths
  Umat_Zgamma = matrix(0, nr = nresp*nscns*nalts, nc=ndraw)

    
  #loop over respondents
  for(resp in 1:nresp){
    #pull betas
    betadraw <- hmnl_fit$Beta
    #transpose dimensions of betadraw array
    betadraw <- aperm(betadraw, c(2,1,3))
    #convert to matrix for faster mean over respondents
    betadraw <- matrix(betadraw, nr= dim(betadraw)[1])
    #take mean over respondents to get mean posterior betas for this model
    post_betabar_mat<- matrix(colMeans(betadraw), nr=ndraw)
      
    #use upper level to get mean of dist of heterogeneity
    gammadraw=hmnl_fit$Gamma
    #transpose dimensions of gammadraw array
    gammadraw <- aperm(gammadraw, c(2,1,3))
    #multiply by Z to get mean of dist of het
    Zgamma_mat <- matrix(Z[resp,]%*%gammadraw, nr=ndraw)

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
  #save probabilities
  probs_postbetas <- (Umat_postbetas/sums)
    
   #find probabilities for each task, resp and draw for Zgamma
  Umat_Zgamma <- matrix(Umat_Zgamma, nr = nalts) 
  sums <- t(matrix(rep(colSums(Umat_Zgamma),nalts), nc=nalts))
  #save probabilities
  probs_Zgamma <- (Umat_Zgamma/sums)
  
  #find location of highest prob postbetas
  locs_postbetas <- apply(probs_postbetas,2,which.max)
  
  #find location of highest prob Zgamma
  locs_Zgamma <- apply(probs_Zgamma,2,which.max)
  
  #calculate hits postbetas
  hits_postbetas <- double(nresp*nscns*ndraw)
  hits_postbetas[locs_postbetas==rep(test_Y_stacked,ndraw)] <- 1
  
  #calculate hit probs postbetas
  hit_probs_postbetas <- colSums(probs_postbetas*diag(nalts)[,rep(test_Y_stacked, ndraw)])
  
  #calculate hits Zgamma
  hits_Zgamma <- double(nresp*nscns*ndraw)
  hits_Zgamma[locs_Zgamma==rep(test_Y_stacked, ndraw)] <- 1
  
  #calculate hit probs Zgamma
  hit_probs_Zgamma <- colSums(probs_Zgamma*diag(nalts)[,rep(test_Y_stacked, ndraw)])
  
  hitrates=c(mean(hits_postbetas),mean(hits_Zgamma))
  hitprobs=c(mean(hit_probs_postbetas), mean(hit_probs_Zgamma))
  
  
  return(list(hit_prob=hitprobs, hit_rate=hitrates))
}
