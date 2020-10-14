# Start with loo.

hit_rate = function(betadraw,Y,X) {
  # Compute the hit rate for the indicated model.
  #   betadraw - post-burn-in beta draws for the model of interest
  #   Y - choices (hold-out tasks or sample)
  #   X - design matrices (hold-out tasks or sample)
  
  ndraw = length(betadraw[1,1,])         # Number of draws.
  nresp = length(Y)                      # Number of respondents.
  nscns = length(Y[[1]])                 # Number of hold-out choice tasks.
  nalts = length(X[[1]][,1])/nscns       # Number of alternatives in each hold-out choice task.
  
  # RESPONDENT LOOP: Average the predicted probability across respondents.
  hit_list = NULL
  for (resp in 1:nresp) {
    
    # HOLD-OUT LOOP: Average the hits across hold-out scenarios.
    for (scn in 1:nscns) {
      Y_scn = Y[[resp]][scn]
      X_scn = X[[resp]][((nalts*scn)-(nalts-1)):(nalts*scn),]
      
      # MCMC-DRAW LOOP: Average the hits across post-burn-in MCMC draws.
      Y_predict = NULL
      for (draw in 1:ndraw) {
        
        # Predict Y.
        Y_predict = rbind(Y_predict,which.max(X_scn%*%betadraw[resp,,draw]))
        
      }
      
      # Compare the prediction with the hold-out response for this scenario.
      hit_list = cbind(hit_list,(rep(Y_scn,ndraw)==as.vector(Y_predict))*1)
      
    }
    print(resp)
  }
  return(round(sum(hit_list)/length(hit_list),3))
}

hit_prob = function(betadraw,Y,X) {
  # Compute the hit probability for the indicated model.
  #   betadraw - post-burn-in beta draws for the model of interest
  #   Y - choices (hold-out tasks or sample)
  #   X - design matrices (hold-out tasks or sample)
  
  ndraw = length(betadraw[1,1,])         # Number of draws.
  nresp = length(Y)                      # Number of respondents.
  nscns = length(Y[[1]])                 # Number of hold-out choice tasks.
  nalts = length(X[[1]][,1])/nscns       # Number of alternatives in each hold-out choice task.
  
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
  
  # RESPONDENT LOOP: Average the predicted probability across respondents.
  resp_prob = matrix(NA,nrow=nresp,ncol=1); test = NULL
  for (resp in 1:nresp) {
    
    # HOLD-OUT LOOP: Average the predicted probability across hold-out scenarios.
    hold_out_prob = matrix(NA,nrow=nscns,ncol=1); test2 = NULL
    for (scn in 1:nscns) {
      Y_scn = Y[[resp]][scn]
      X_scn = X[[resp]][((nalts*scn)-(nalts-1)):(nalts*scn),]
      
      # MCMC-DRAW LOOP: Average the predicted probability across post-burn-in MCMC draws.
      mcmc_draw_prob = matrix(NA,nrow=ndraw,ncol=1)
      for (draw in 1:ndraw) {
        
        # Compute the likelihood of making the observed choice for the given scenario using the given model.
        mcmc_draw_prob[draw] = exp(ll_mnl(betadraw[resp,,draw],Y_scn,X_scn))
        
      }
      hold_out_prob[scn] = sum(mcmc_draw_prob)/ndraw
      test2 = cbind(test2,sum(mcmc_draw_prob)/ndraw)
      
    }
    resp_prob[resp,] = sum(hold_out_prob)/nscns
    test = cbind(test,test2)
    print(c(resp,resp_prob[resp]),digits=2)
  }
  return(round(sum(resp_prob)/nresp,3))
}
