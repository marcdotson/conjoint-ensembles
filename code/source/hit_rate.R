hit_rate = function(betadraw,Y,X) {
  # Compute the hit rate for the indicated model.
  #   betadraw - post-burn-in beta draws for the model of interest
  #   Y - choices (hold-out tasks or sample)
  #   X - design matrices (hold-out tasks or sample)
  
  ndraw = length(betadraw[,1,1])         # Number of draws.
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
