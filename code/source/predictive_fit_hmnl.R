predictive_fit_hmnl = function(hmnl_draws, test_X, test_Y, test_Z){
  # Compute the hit rate for the indicated model.
  #   hmnl_fit - hmnl output with log_lik, betadraws, gammadraws, and Omegadraws
  #   test_Y - choices (hold-out sample)
  #   test_X - design matrices (hold-out sample)
  #   test_Z - matrix of covariates

  ndraw <- length(hmnl_draws$Gamma[,1,1]) # Number of draws
  nresp <- length(test_Y[,1])           # Number of respondents
  nscns <- length(test_X[1, ,1,1])      # Number of choice tasks
  nalts <- length(test_X[1,1, ,1])      # Number of alternatives 
  nlvls <- length(test_X[1,1,1, ])      # Number of att levels
  if(is.null(test_Z)){test_Z <- matrix(1, nr=nresp, nc = 1)}
  
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
    #using mean (over draws) of the mean of the post dist
    Umat <-  matrix(0, nr = nresp*nscns*nalts)
    
  #loop over respondents
  for(resp in 1:nresp){
    #use upper level to get mean of dist of heterogeneity
    gammadraws=hmnl_draws$Gamma
    #transpose dimensions of gammadraw array
    meangammas <- apply(gammadraws,c(2,3),mean)
    #multiply by Z to get mean of dist of het
    betas <- matrix(test_Z[resp,]%*%meangammas, nc=1)

    #get utility for each alternative
    # Umat_meangammas[((resp-1)*nalts*nscns+1):((resp)*nalts*nscns),] <- 
    #   exp(test_X_stacked[((resp-1)*nalts*nscns+1):((resp)*nalts*nscns),]%*%
    #         matrix(betas))
    Umat[((resp-1)*nalts*nscns+1):((resp)*nalts*nscns),] <- 
      exp(test_X_stacked[((resp-1)*nalts*nscns+1):((resp)*nalts*nscns),]%*%
            matrix(betas))
  }
    
  #find probabilities for each task, resp
  # Umat_byscn <- matrix(Umat_meangammas, nr = nalts) 
  Umat_byscn <- matrix(Umat, nr = nalts) 
  sums <- t(matrix(rep(colSums(Umat_byscn),nalts), nc=nalts))
  probs <- (Umat_byscn/sums)
  
  #find location of highest prob
  locs <- apply(probs,2,which.max)
  
  #calculate hits meangammas
  hits <- double(nresp*nscns*ndraw)
  hits[locs ==test_Y_stacked] <- 1
  
  #calculate hit probs meangammas
  hit_probs<- colSums(probs*diag(nalts)[,test_Y_stacked])
  
  return(list(hit_prob=mean(hit_probs), hit_rate=mean(hits)))
}
