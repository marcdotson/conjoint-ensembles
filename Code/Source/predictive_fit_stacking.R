predictive_fit_stacking <- function(member_draws, validate_X, validate_Z){
  # Compute predictions for the stacking meta-learner for the given ensemble member.
  #   member_draws - hmnl output with Gamma and Sigma draws
  #   validate_X - design matrices (hold-out validation sample)
  #   validate_Z - matrix of covariates (hold-out validation sample)
  
  nresp <- dim(validate_X)[1]      # Number of respondents
  nscns <- dim(validate_X)[2]      # Number of choice tasks
  nalts <- dim(validate_X)[3]      # Number of alternatives 
  nlvls <- dim(validate_X)[4]      # Number of att levels
  if(is.null(validate_Z)) { validate_Z <- matrix(1, nr=nresp, nc = 1) }
  
  #stack resps and scns to avoid loops (this needs changed if using hold out tasks)
  validate_X_stacked <- NULL
  for(resp in 1:nresp) {
    for(scn in 1:nscns) {
      validate_X_stacked <- rbind(validate_X_stacked, validate_X[resp, scn,,])
    }
  }
  
  #get utilities for computing predictions for the ensemble member
  #using mean of the post dist (already computed during ensemble fit)
  Umat <-  matrix(0, nr = nresp * nscns * nalts)
  
  #loop over respondents
  for(resp in 1:nresp) {
    #use upper level draws as mean of dist of heterogeneity
    gammadraws <- member_draws$Gamma
    
    #multiply by Z to get mean of dist of het
    betas <- matrix(validate_Z[resp,] %*% gammadraws, nc = 1)
    
    #get utility for each alternative
    Umat[((resp - 1) * nalts * nscns + 1):((resp) * nalts * nscns),] <- 
      exp(
        validate_X_stacked[((resp - 1) * nalts * nscns + 1):((resp) * nalts * nscns),] %*% matrix(betas)
      )
  }
  
  #find probabilities for each task, resp
  Umat_byscn <- matrix(Umat, nr = nalts)
  sums <- t(matrix(rep(colSums(Umat_byscn), nalts), nc = nalts))
  probs <- (Umat_byscn / sums)
  
  #find location of highest prob
  locs <- apply(probs, 2, which.max)
  predicted_Y <- matrix(locs, nrow = nresp, ncol = nscns, byrow = TRUE)
  
  return(predicted_Y)
}
