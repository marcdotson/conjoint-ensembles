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

model_fit = function(betadraw,Y,X) {
  # Compute the LMD, hit rates, and hit probability for the indicated model.
  #   betadraw - post-burn-in beta draws for the model of interest
  #   Y - choices (hold-out tasks or sample)
  #   X - design matrices (hold-out tasks or sample)
  
  ndraw = length(betadraw[1,1,])         # Number of draws.
  nresp = length(Y)                      # Number of respondents.
  nscns = length(Y[[1]])                 # Number of hold-out choice tasks.
  nalts = length(X[[1]][,1])/nscns       # Number of alternatives in each hold-out choice task.
  
  # RESPONDENT LOOP: Compute the LL, predict Y, and average the predicted probability across respondents.
  llike = array(NA,dim=c(ndraw,nscns,nresp)); hits = NULL; resp_prob = matrix(NA,nrow=nresp,ncol=1)
  for (resp in 1:nresp) {
    
    # HOLD-OUT LOOP: Compute the LL, predict Y, and average the predicted probability across hold-out scenarios.
    hold_out_prob = matrix(NA,nrow=nscns,ncol=1)
    for (scn in 1:nscns) {
      Y_scn = Y[[resp]][scn]
      X_scn = X[[resp]][((nalts*scn)-(nalts-1)):(nalts*scn),]
      
      # MCMC-DRAW LOOP: Compute the LL, predict Y, and average the predicted probability across post-burn-in MCMC draws.
      Y_predict = NULL; mcmc_draw_prob = matrix(NA,nrow=ndraw,ncol=1)
      for (draw in 1:ndraw) {
        mcmc_draw_llike = ll_mnl(betadraw[resp,,draw],Y_scn,X_scn)           # Compute the LL.
        llike[draw,scn,resp] = mcmc_draw_llike                               # Store the LL.
        Y_predict = rbind(Y_predict,which.max(X_scn%*%betadraw[resp,,draw])) # Store the predicted Y.
        mcmc_draw_prob[draw,] = exp(mcmc_draw_llike)                         # Store the predicted probability.
      }
      hits = cbind(hits,(rep(Y_scn,ndraw)==as.vector(Y_predict))*1) # Identify hits.
      hold_out_prob[scn] = sum(mcmc_draw_prob)/ndraw                # Average the predicted probabilities.
    }
    resp_prob[resp,] = sum(hold_out_prob)/nscns          # Average the predicted probabilities.
    print(paste(round(resp/nresp,2)*100,"Percent Done")) # Print progress.
  }
  print(paste("Hit Rate:",round(sum(hits)/length(hits),3),"Hit Prob:",round(sum(resp_prob)/nresp,3)))
  return(list(round(sum(hits)/length(hits),3),round(sum(resp_prob)/nresp,3),llike))
}

