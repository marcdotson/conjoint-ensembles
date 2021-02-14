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
  if(is.null(Z)){Z <- matrix(1, nr=nresp, nc = 1)}
  
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
    #1) using mean of the distribution of heterogeneity as predicted part worths
    Umat_gammadraws = matrix(0, nr = nresp*nscns*nalts, nc=ndraw)
    #2) using mean (over draws) of the mean of the post dist
    Umat_meangammas = matrix(0, nr = nresp*nscns*nalts)
    
  #loop over respondents
  for(resp in 1:nresp){
    #use upper level to get mean of dist of heterogeneity
    gammadraws=hmnl_fit$Gamma
    #transpose dimensions of gammadraw array
    gammadraws <- aperm(gammadraws, c(2,1,3))
    #multiply by Z to get mean of dist of het
    gammadraws_mat <- matrix(Z[resp,]%*%gammadraws, nr=ndraw)

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
  probs_gammadraws <- (Umat_gammadraws_byscn/sums)
  #find probabilities for each task, resp and draw for meangammas
  Umat_meangammas_byscn <- matrix(Umat_meangammas, nr = nalts) 
  sums <- t(matrix(rep(colSums(Umat_meangammas_byscn),nalts), nc=nalts))
  probs_meangammas <- (Umat_meangammas_byscn/sums)
  
  #find location of highest prob gammadraws
  locs_gammadraws <- apply(probs_gammadraws,2,which.max)
  #find location of highest prob meangammas
  locs_meangammas <- apply(probs_meangammas,2,which.max)
  
  #calculate hits gammadraws
  hits_gammadraws <- double(nresp*nscns*ndraw)
  hits_gammadraws[locs_gammadraws==rep(test_Y_stacked, ndraw)] <- 1
  #calculate hits meangammas
  hits_meangammas <- double(nresp*nscns*ndraw)
  hits_meangammas[locs_meangammas==test_Y_stacked] <- 1
  
  #calculate hit probs gammadraws
  hit_probs_gammadraws<- colSums(probs_gammadraws*diag(nalts)[,rep(test_Y_stacked, ndraw)])
  #calculate hit probs meangammas
  hit_probs_meangammas<- colSums(probs_meangammas*diag(nalts)[,test_Y_stacked])
  
  return(list(hit_prob_gammadraws=mean(hit_probs_gammadraws), hit_rate_gammadraws=mean(hits_gammadraws),
              hit_prob_meangammas=mean(hit_probs_meangammas), hit_rate_meangammas=mean(hits_meangammas)))
}
