
# Load packages.
library(tidyverse)
library(rstan)
library(bayesplot)
library(tidybayes)
library(loo)

# Functions for model fit/comparison

predictive_fit_ensemble = function(ensemble_weights, ensemble_fit, test_x, test_y){
  # Compute the hit rate, hit prob, and loo metrics for the ensemble model.
  #   ensemble_weights - estimated weights for each of the models
  #   ensemble_fit - ensemble output with log_lik and betadraws for each model
  #   test_Y - choices (hold-out sample)
  #   test_X - design matrices (hold-out sample)
 
  n_ens <- length(ensemble_fit)
  ndraw <- length(ensemble_fit[[1]]$Beta[,1,1])         # Number of draws
  nresp <- length(test_Y[,1])           # Number of respondents
  nscns <- length(test_X[1, ,1,1])      # Number of choice tasks
  nalts <- length(test_X[1,1, ,1])      # Number of alternatives 
  nlvls <- length(test_X[1,1,1, ])      # Number of att levels
  
  #weight log_lik for each model to get log_lik for ensemble
  LLarray_ens = array(0,dims)
  for(k in 1:n_ens){
    #extract log_lik array from each stanfit object
    LLarray_ens <- LLarray_ens + weights[k]*ensemble_fit[[k]]$log_lik
  }  
  
  #get effective sample size
  r_eff_ens <- loo::relative_eff(x = exp(LLarray_ens), cores = cores)
  
  #apply PSIS via loo to ensemble likelihoods (loo fit metrics)
  loo_fit_ens <- loo::loo.array(LLarray, r_eff = r_eff,
                                cores = cores, save_psis = FALSE)
  
  
  
  #stack resps and scns to avoid loops (this needs changed if using hold out tasks)
  test_X_stacked <- NULL
  for(resp in 1:nresp){
    for(scn in 1:nscns){
      test_X_stacked <- rbind(test_X_stacked,test_X[resp,scn,,])
    }
  }
  
  #stack scn choices to avoid loops
  test_Y_stacked <- matrix(t(test_Y),nr=1)
  
  #loop over ensemble models to calculate individual hit rates
  hit_rate_vec=double(n_ens)
  hit_prob_vec=double(n_ens)
  for(n in 1:n_ens){
    betadraw <- ensemble_fit[[n]]$Beta
    
    #average betadraws over in-sample respondents for hold-out sample (no upper level)
    betadraw_avg <- matrix(double(ndraw*nlvls), nc=nlvls)
    for(draw in 1:ndraw){
      betadraw_avg[draw,] <- colMeans(betadraw[draw,,])
    }
    
    #find probabilities for each
    Umat <- exp(test_X_stacked%*%t(betadraw_avg))
    Umat <- matrix(Umat, nr = 3) 
    sums <- t(matrix(rep(colSums(Umat),3), nc=3))
    probs <- Umat/sums
    
    #find location of highest prob
    locs <- apply(probs,2,which.max)
    
    #identify hitrate and store
    hits <- double(length(locs))
    hits[ locs== rep(test_Y_stacked,ndraw) ] <- 1
    hit_rate_vec[n]=mean(hits)
    
    #identify hitprobs and store
    prob_select <- colSums(probs*diag(nalts)[,rep(test_Y_stacked,ndraw)])
    hit_prob_vec[n]=mean(prob_select)
  }
  
  hit_rate_ens <- hit_rate_vec%*%ensemble_weights
  hit_prob_ens <- hit_prob_vec%*%ensemble_weights
  return(list(hit_prob=hit_prob_ens, hit_rate=hit_rate_ens, loo_fit=loo_fit_ens))
}

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
