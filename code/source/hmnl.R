hier_mnl_e = function (Data, Prior, Mcmc, Cont) {
  # SUMMARY: Implements a Random-Walk Metropolis for a hierarchical
  # MNL with a multivariate normal distribution of heterogeneity.
  
  # FUNCTION ARGUMENTS:
  # Data - list(nalts,data,Z)
  #  nalts - number of alternatives
  #  data - list of lists with data for each respondent
  #   data$yy - vector of multinomial outcomes
  #   data$XX - design matrix
  #  Z - covariates for the upper-level model
  # Prior - list(gammabar,Agamma,nu,V)
  #  gammabar - means for normal prior on Gamma
  #  Agamma - precision matrix for normal prior on Gamma
  #  nu - df for IW on prior on Vbeta
  #  V - location for IW prior on Vbeta
  # Mcmc - list(R,keep,step)
  #  R - number of draws
  #  keep - thinning parameter
  #  step - RW step (scaling factor) for the beta draws
  #  sim_ind - indicates a simulation experiment
  #  cont_ind - indicates a run continuation
  # Cont - list(out_bstep,out_oldbetas,out_oldgamma,out_oldVbeta)
  #  out_bstep - optional ending step from previous run
  #  out_oldbetas - optional oldbetas from previous run
  #  out_oldgamma - optional oldgamma from previous run
  #  out_oldVbeta - optional oldVbeta from previous run
  
  # FUNCTION OUTPUT:
  # betadraw - posterior draws of individual-level betas
  # gammadraw - posterior draws of the upper level coefficient matrix
  # Vbetadraw - posterior draws of the upper level covariance matrix
  # llikedraw - draws of the log likelihood
  # baccept - acceptance rate for the beta step
  # bstepkeep - final adjust beta step after burn-in
  # Cont - list(out_bstep,out_oldbetas,out_oldgamma,out_oldVbeta)
  #  out_bstep - ending step to continue from
  #  out_oldbetas - oldbetas to continue from
  #  out_oldgamma - oldgamma to continue from
  #  out_oldVbeta - oldVbeta to continue from
  
  # Assign values from the function arguments.
  nalts = Data$nalts
  data = Data$data
  nresp = length(data)
  Z = Data$Z
  nz = ncol(Data$Z)
  gammabar = Prior$gammabar
  nvars = ncol(Prior$gammabar)
  Agamma = Prior$Agamma
  nu = Prior$nu
  V = Prior$V
  R = Mcmc$R
  keep = Mcmc$keep
  print = Mcmc$print
  sim_ind = Mcmc$sim_ind
  cont_ind = Mcmc$cont_ind
  
  # Define matrices for posterior draws, set initial clock time.
  betadraw = array(double(floor(R/keep)*nresp*nvars),dim=c(nresp,nvars,floor(R/keep)))
  gammadraw = matrix(double(floor(R/keep)*nvars*nz),ncol=nvars*nz)
  Vbetadraw = matrix(double(floor(R/keep)*nvars*nvars),ncol=nvars*nvars)
  llikedraw = array(0,dim=c(floor(R/keep), nresp)) 
  baccept = array(0,dim=c(R/keep))
  bstepkeep = array(0,dim=c(R/keep))
  itime = proc.time()[3]
  cat("MCMC Iteration (estimated time to end in hours | step | baccept | ll )",fill=TRUE)
  
  # Setup for RW chain.
  if (cont_ind == 0) {
    step = Mcmc$step
    if (sim_ind==0) oldbetas = matrix(double(nresp*nvars),ncol=nvars)
    if (sim_ind==1) oldbetas = matrix(Data$beta_true,ncol=nvars)
    oldgamma = matrix(double(nvars*nz),ncol=nvars)
    oldVbeta = diag(nvars)
    oldVbetai = diag(nvars)
  }
  if (cont_ind == 1) {
    step = Cont$out_bstep
    oldbetas = Cont$out_oldbetas
    oldgamma = matrix(Cont$out_oldgamma,ncol=nvars)
    oldVbeta = Cont$out_oldVbeta
    oldVbetai = Cont$out_oldVbeta
  }
  
  # Run the RW chain to generate as many draws from the posterior as specified in R.
  for (rep in 1:R) {
    # Initial log likelihood values for each iteration.
    logold = lognew = 0
    loglike = 0
    
    bnaccept = 0
    for (resp in 1:nresp) {
      # Beta old and candidate draws.
      betad = oldbetas[resp,]
      betac = as.vector(rmvnorm(1,mean=betad,sigma=(step*oldVbeta)))
      
      # Log likelihood with the old gamma draws and old/candidate beta draws.
      logold = llmnl(betad,data[[resp]]$yy,data[[resp]]$XX)
      lognew = llmnl(betac,data[[resp]]$yy,data[[resp]]$XX)
      
      # Log of the MVN distribution of heterogeneity over all betas.
      loghold = -0.5*(t(betad)-Z[resp,]%*%oldgamma)%*%oldVbetai%*%(betad-t(Z[resp,]%*%oldgamma))
      loghnew = -0.5*(t(betac)-Z[resp,]%*%oldgamma)%*%oldVbetai%*%(betac-t(Z[resp,]%*%oldgamma))
      
      # Beta log posteriors.
      lpostold = logold + loghold
      lpostnew = lognew + loghnew
      
      # Compare the old and candidate posteriors and compute alpha (second-stage prior cancels out).
      diff = exp((lpostnew) - (lpostold))
      if (diff == "NaN" || diff == Inf) {
        alpha = -1 # If the number doesn't exist, always reject.
      } else {
        alpha = min(1,diff)
      }
      unif = runif(1)
      if (unif < alpha) {
        oldbetas[resp,] = betac
        bnaccept = bnaccept + 1
        loglike = loglike + lognew
        if(rep%%keep==0) llikedraw[rep/keep, resp] = lognew
      } else {
        loglike = loglike + logold
        if(rep%%keep==0) llikedraw[rep/keep, resp] = logold
      }
    }
    
    # Gamma and Vbeta draw (distribution of heterogeneity over beta).
    out = rmultireg(oldbetas,Z,gammabar,Agamma,nu,V)
    oldgamma = out$B
    oldVbeta = out$Sigma
    oldVbetai = chol2inv(chol(oldVbeta))
    
    # Modify the RW step sizes to constrain acceptance rates during burn-in (R/3).
    if (rep%%5 == 0 & cont_ind == 0) {
      if (rep < (R/2)) {
        # Update step.
        if (bnaccept/nresp < .20) {
          step = step*0.95
        }
        if (bnaccept/nresp > .60) {
          step = step*1.05
        }
      }
    }
    
    # Print progress.
    if (rep%%(keep*5) == 0) {
      ctime = proc.time()[3]
      timetoend = ((ctime - itime)/rep)*(R - rep)
      bacceptr=bnaccept/nresp
      cat(" ",rep," (",round((timetoend/60)/60,2),"|",round(step,5),"|",round(bacceptr,2),"|",round(loglike,2),")",fill = TRUE)
    }
    
    # Print chart less often.
    if (rep%%print == 0) {
      par(mfrow=c(2,1))
      if (sim_ind==0) { plot(rowSums(llikedraw),type="l"); matplot(gammadraw,type="l") }
      if (sim_ind==1) {
        plot(rowSums(llikedraw),type="l")
        matplot(gammadraw[,1:nz],type="l",col=c(1:nz)); abline(h=Data$Gamma[,1],col=c(1:nz))
      }
    }
    
    # Save the posterior draws.
    mkeep = rep/keep
    if (mkeep*keep == (floor(mkeep)*keep)) {
      betadraw[,,mkeep] = oldbetas
      gammadraw[mkeep,] = as.vector(oldgamma)
      Vbetadraw[mkeep,] = as.vector(oldVbeta)
      baccept[mkeep] = bnaccept/nresp
      bstepkeep[mkeep] = step
    }
    
    # Save out continuation files.
    if (rep%%R == 0) {
      Cont = list(out_oldbetas = betadraw[,,R/keep],out_oldgamma = matrix(gammadraw[R/keep,],byrow=TRUE,ncol=(nvars*nz)),
                  out_oldVbeta = matrix(Vbetadraw[R/keep,],byrow=TRUE,ncol=nvars),out_bstep = step)
    }
  }
  
  # Print total run time.
  ctime = proc.time()[3]
  cat(" Total Time Elapsed (in Hours): ",round(((ctime - itime)/60)/60,2),fill = TRUE)
  
  # Output.
  return(list(betadraw=betadraw,gammadraw=gammadraw,Vbetadraw=Vbetadraw,
              llikedraw=llikedraw,baccept=baccept,bstepkeep=bstepkeep,Cont=Cont))
}
