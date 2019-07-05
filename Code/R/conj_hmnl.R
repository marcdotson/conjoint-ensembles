conj_hmnl = function (Data, Prior, Mcmc, Cont) {
  # This function implements a random-walk Metropolis-Hastings algorithm for a hierarchical
  # MNL with a multivariate normal distribution of heterogeneity and a conjuctive screen
  # based on product relevance as determined by attribute perceptions.
  
  # Describe and Assign Function Arguments ----------------------------------
  # Data = list(data,W,Q,D).
  data = Data$data                                        # List of lists with y, X, S.
  W = Data$W                                              # Covariates for the upper-level preference model.
  if (Mcmc$het_ind==1) {
    Q = Data$Q                                            # Covariates for the upper-level screening model.
    D = Data$D                                            # "Design matrix" for the upper-level screening model.
  }
  if (Mcmc$sim_ind==1) {
    Beta = Data$Beta                                      # True values of Beta for parameter recovery.
    Tau = Data$Tau                                        # True values of Tau for parameter recovery.
    if (Mcmc$het_ind==1) Theta = Data$Theta               # True values of Theta for parameter recovery.
    Gamma = Data$Gamma                                    # True values of Gamma for parameter recovery.
    Vbeta = Data$Vbeta                                    # True values of Vbeta for parameter recovery.
    if (Mcmc$het_ind==0) theta = Data$theta               # True values of theta for parameter recovery.
    if (Mcmc$het_ind==1) Omega = Data$Omega               # True values of Omega for parameter recovery.
  }
  
  # Prior = list(gammabar,Agamma,nu,V,omegabar,Vomegabar).
  gammabar = Prior$gammabar                               # Means for normal prior on Gamma.
  Agamma = Prior$Agamma                                   # Precision matrix for normal prior on Gamma.
  nu = Prior$nu                                           # DF for IW prior on Vbeta.
  V = Prior$V                                             # Location for IW prior on Vbeta.
  if (Mcmc$het_ind==1) {
    omegabar = Prior$omegabar                             # Means for normal prior on Omega.
    Vomegabar = Prior$Vomegabar                           # Covariance matrix for normal prior on Omega.
  }
  
  # Mcmc = list(R,keep,bstep,ostep,sim_ind,out_ind,het_ind,cont_ind).
  R = Mcmc$R                                              # Number of iterations in the Markov chain.
  keep = Mcmc$keep                                        # Thinning parameter.
  bstep = Mcmc$bstep                                      # RW step (scaling factor) for the beta draws.
  if (Mcmc$het_ind==1) ostep = Mcmc$ostep                 # RW step (scaling factor) for the Omega draws.
  sim_ind = Mcmc$sim_ind                                  # Indicates a simulation experiment.
  out_ind = Mcmc$out_ind                                  # Indicates an outside good.
  het_ind = Mcmc$het_ind                                  # Indicates heterogeneity in Theta.
  cont_ind = Mcmc$cont_ind                                # Indicates a run continuation.
  
  # Assign values from the function arguments.
  nresp = length(data)                                    # Number of respondents.
  nscns = length(data[[1]]$y)                             # Number of choice tasks.
  nalts = length(data[[1]]$X[,1])/nscns                   # Number of alternatives in each choice task.
  nlvls = ncol(data[[1]]$S)                               # Number of attribute levels (including the brand).
  nvars = ncol(data[[1]]$X)                               # Number of attribute levels (excluding reference levels).
  npcov = ncol(W)                                         # Number of covariates for preference heterogeneity.
  if (het_ind==1) nscov = ncol(Q)-1                       # Number of covariates for screening heterogeneity.
  
  # Describe and Initialize Function Output ---------------------------------
  # Respondent-level parameter draws.
  betadraw = array(double(floor(R/keep)*nresp*nvars),dim=c(nresp,nvars,floor(R/keep)))
  taudraw = array(double(floor(R/keep)*nresp*nlvls),dim=c(nresp,nlvls,floor(R/keep)))
  if (het_ind==1) thetadraw = array(double(floor(R/keep)*nresp*nlvls),dim=c(nresp,nlvls,floor(R/keep)))
  
  # Aggregate-level parameter draws.
  Gammadraw = matrix(double(floor(R/keep)*nvars*npcov),ncol=nvars*npcov)
  Vbetadraw = matrix(double(floor(R/keep)*nvars*nvars),ncol=nvars*nvars)
  if (het_ind==0) thetadraw = matrix(double(floor(R/keep)*nlvls),ncol=nlvls)
  if (het_ind==1) Omegadraw = matrix(double(floor(R/keep)*(nscov+1)*nlvls),ncol=(nscov+1)*nlvls)
  
  # Diagnostic draws and initial clock time.
  llikedraw = double(floor(R/keep))                        # Log likelihood.
  baccept = array(0,dim=c(R/keep))                         # Beta acceptance rate.
  if (het_ind==1) oaccept = array(0,dim=c(R/keep))         # Omega acceptance rate.
  bstepdraw = array(0,dim=c(R/keep))                       # RW step adjusted during burn-in.
  if (het_ind==1) ostepdraw = array(0,dim=c(R/keep))       # RW step adjusted during burn-in.
  itime = proc.time()[3]                                   # Initial clock time.
  
  # Initialize MCMC ---------------------------------------------------------
  if (het_ind==0) cat("MCMC Iteration (estimated time to end in hrs/min | bstep | baccept | ll )",fill=TRUE)
  if (het_ind==1) cat("MCMC Iteration (estimated time to end in hrs/min | bstep | baccept | ostep | oaccept | ll )",fill=TRUE)
  
  # Initialize values.
  if (cont_ind == 0) {
    if (sim_ind==0) {
      oldbetas = matrix(double(nresp*nvars),ncol=nvars)
      oldtaus = matrix(double(nresp*nlvls),ncol=nlvls)
    }
    if (sim_ind==1) {
      oldbetas = matrix(Beta,ncol=nvars)
      oldtaus = matrix(Tau,ncol=nlvls)
    }
    if (het_ind==1) {
      if (sim_ind==0) oldtheta = matrix(rep(.5,nresp*nlvls),ncol=nlvls)
      if (sim_ind==1) oldtheta = matrix(Theta,ncol=nlvls)
    }
    canSigma = 0.1 * diag(nvars)
    oldGamma = matrix(double(nvars*npcov),ncol=nvars)
    oldVbeta = diag(nvars)
    oldVbetai = diag(nvars)
    if (het_ind==0) oldtheta = matrix(rep(.5,nlvls),ncol=nlvls)
    if (het_ind==1) {
      oldOmega = matrix(double(nlvls*(nscov+1)),ncol=nscov+1)
      Vomegabari = chol2inv(chol(Vomegabar))
    }
  }
  
  # Initialize values and use the previous draws for continued runs.
  if (cont_ind == 1) {
    bstep = Cont$out_bstep
    if (het_ind==1) ostep = Cont$out_ostep
    oldbetas = Cont$out_oldbetas
    oldtaus = Cont$out_oldtaus
    if (het_ind==1) oldtheta = matrix(Cont$out_oldtheta,ncol=nlvls)
    canSigma = 0.1 * diag(nvars)
    oldGamma = matrix(Cont$out_oldGamma,ncol=nvars)
    oldVbeta = Cont$out_oldVbeta
    oldVbetai = Cont$out_oldVbeta
    if (het_ind==0) oldtheta = Cont$out_oldtheta
    if (het_ind==1) {
      oldOmega = matrix(Cont$out_oldOmega,ncol=nscov+1)
      Vomegabari = chol2inv(chol(Vomegabar))
    }
  }
  
  # Run the MCMC ------------------------------------------------------------
  # The Markov chain will run for R iterations.
  sourceCpp("Source/conj_hmnl.cpp")
  for (rep in 1:R) {
    # Initial log likelihood values for each iteration.
    logold = lognew = 0
    loglike = 0
    
    # Respondent-level loop.
    bnaccept = 0; if (het_ind==1) onaccept = 0
    for (resp in 1:nresp) {
      # Beta old and candidate draws.
      betad = oldbetas[resp,]
      betac = as.vector(rmvnorm(1,mean=betad,sigma=(bstep*canSigma)))
      # if (rep < (R/3)) betac = as.vector(rmvnorm(1,mean=betad,sigma=(bstep*oldVbeta)))
      # if (rep >= (R/3)) betac = as.vector(rmvnorm(1,mean=betad,sigma=(bstep*Vbeta_fixed)))

      # Evaluate the likelihood for old and candidate beta draws.
      out = conj_llmnl_Cpp(nscns,nalts,nlvls,nvars,as.vector(data[[resp]]$y),data[[resp]]$X,data[[resp]]$X,
                           data[[resp]]$S,oldtaus[resp,],betad,betac,out_ind)
      logold <- out[[1]]
      lognew <- out[[2]]

      # Log of the MVN distribution of heterogeneity over all betas.
      loghold = -0.5*(t(betad)-W[resp,]%*%oldGamma)%*%oldVbetai%*%(betad-t(W[resp,]%*%oldGamma))
      loghnew = -0.5*(t(betac)-W[resp,]%*%oldGamma)%*%oldVbetai%*%(betac-t(W[resp,]%*%oldGamma))

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
      } else {
        loglike = loglike + logold
      }
      
      # Draw tau.
      oldtaus[resp,] = draw_tau_Cpp(resp,nresp,nscns,nalts,nlvls,as.vector(data[[resp]]$y),data[[resp]]$S,
                                    data[[resp]]$X%*%oldbetas[resp,],oldtaus[resp,],oldtheta,out_ind,het_ind)
    }
    
    # Draw Gamma and Vbeta (distribution of heterogeneity over beta).
    out = rmultireg(oldbetas,W,gammabar,Agamma,nu,V)
    oldGamma <- out$B
    oldVbeta <- out$Sigma
    oldVbetai <- chol2inv(chol(oldVbeta))
    
    # # Fix Vbeta for post-burn-in beta proposal density.
    # if (rep == floor(R/3)) Vbeta_fixed = oldVbeta
    
    if (het_ind==0) {
      # Draw theta.
      oldtheta = matrix(rbeta(sum(nlvls),colSums(oldtaus)+1,nresp-colSums(oldtaus)+2),nrow=1,ncol=sum(nlvls))
    }
    if (het_ind==1) {
      # Draw Omega and theta.
      out = draw_omega_theta_Cpp(nresp,nlvls,nscov,Q,D,oldOmega,ostep,onaccept,omegabar,Vomegabar,Vomegabari,ifelse(oldtaus==0,1,0),oldtheta)
      oldOmega <- out[[1]]
      onaccept <- out[[2]]
      oldtheta <- out[[3]]
    }
    
    # Houskeeping and Output --------------------------------------------------
    # Modify the RW step sizes to constrain the acceptance rate during burn-in (R/3).
    if (rep%%5 == 0 & cont_ind == 0) {
      if (rep < (R/2)) {
        # Update bstep.
        if (bnaccept/nresp < .20) {
          bstep = bstep*0.95
        }
        if (bnaccept/nresp > .60) {
          bstep = bstep*1.05
        }
        if (het_ind==1) {
          # Update ostep.
          if (onaccept/nlvls < .20) {
            ostep = ostep*0.95
          }
          if (onaccept/nlvls > .60) {
            ostep = ostep*1.05
          }
        }
      }
    }
    
    # Print progress.
    if (rep%%100 == 0) {
      ctime = proc.time()[3]
      timetoend = ((ctime - itime)/rep)*(R - rep)
      bacceptr = bnaccept/nresp
      if (het_ind==1) oacceptr = onaccept/nlvls
      if (het_ind==0) {
        cat(" ",rep," ( hrs/min to end",round((timetoend/60)/60,2),"/",round(timetoend/60,2),"| bstep",
            round(bstep,5),"| baccept",round(bacceptr,2),"| ll",round(loglike,2),")",fill = TRUE)
      }
      if (het_ind==1) {
        cat(" ",rep," ( hrs/min to end",round((timetoend/60)/60,2),"/",round(timetoend/60,2),"| bstep",
            round(bstep,5),"| baccept",round(bacceptr,2),"| ostep",round(ostep,5),"| oaccept",round(oacceptr,2),
            "| ll",round(loglike,2),")",fill = TRUE)
      }
    }
    
    # Print chart less often.
    if (rep%%1000 == 0) {
      par(mfrow=c(3,1))
      if (sim_ind==0) {
        plot(llikedraw,type="l")
        matplot(Gammadraw,type="l",col=1:(nvars*npcov))
        if (het_ind==0) matplot(thetadraw,type="l",col=1:nlvls)
        if (het_ind==1) matplot(Omegadraw,type="l",col=1:((nscov+1)*nlvls))
      }
      if (sim_ind==1) {
        plot(llikedraw,type="l")
        matplot(Gammadraw,type="l",col=1:length(Gamma)); abline(h=Gamma,col=1:length(Gamma))
        if (het_ind==0) { matplot(thetadraw,type="l",col=1:length(theta)); abline(h=theta,col=1:length(theta)) }
        if (het_ind==1) { matplot(Omegadraw,type="l",col=1:length(Omega)); abline(h=Omega,col=1:length(Omega)) }
      }
    }
    
    # Save the posterior draws.
    mkeep = rep/keep
    if (mkeep*keep == (floor(mkeep)*keep)) {
      betadraw[,,mkeep] = oldbetas
      taudraw[,,mkeep] = oldtaus
      if (het_ind==1) thetadraw[,,mkeep] = oldtheta
      Gammadraw[mkeep,] = as.vector(oldGamma)
      Vbetadraw[mkeep,] = as.vector(oldVbeta)
      if (het_ind==0) thetadraw[mkeep,] = as.vector(oldtheta)
      if (het_ind==1) Omegadraw[mkeep,] = as.vector(oldOmega)
      llikedraw[mkeep] = loglike
      baccept[mkeep] = bnaccept/nresp
      if (het_ind==1) oaccept[mkeep] = onaccept/nlvls
      bstepdraw[mkeep] = bstep
      if (het_ind==1) ostepdraw[mkeep] = ostep
    }
    
    # Save out interim continuation files.
    if ( (rep%%round(R/4) == 0)||(rep%%round(R/4)*2 == 0)||(rep%%round(R/4)*3 == 0)) {
      if (het_ind==0) {
        Cont = list(out_oldbetas = betadraw[,,R/keep],
                    out_oldtaus = taudraw[,,R/keep],
                    out_oldGamma = matrix(Gammadraw[R/keep,],byrow=TRUE,ncol=(nvars*npcov)),
                    out_oldVbeta = matrix(Vbetadraw[R/keep,],byrow=TRUE,ncol=nvars),
                    out_oldtheta = matrix(thetadraw[R/keep,],byrow=TRUE,ncol=nlvls),
                    out_bstep = bstep)
      }
      if (het_ind==1) {
        Cont = list(out_oldbetas = betadraw[,,R/keep],
                    out_oldtaus = taudraw[,,R/keep],
                    out_oldtheta = thetadraw[,,R/keep],
                    out_oldGamma = matrix(Gammadraw[R/keep,],byrow=TRUE,ncol=(nvars*npcov)),
                    out_oldVbeta = matrix(Vbetadraw[R/keep,],byrow=TRUE,ncol=nvars),
                    out_oldOmega = matrix(Omegadraw[R/keep,],byrow=TRUE,ncol=((nscov+1)*nlvls)),
                    out_bstep = bstep,out_ostep = ostep)
      }
      save(Cont,file="Interim_Cont.RData")
    }
    
    # Save out continuation files.
    if (rep%%R == 0) {
      if (het_ind==0) {
        Cont = list(out_oldbetas = betadraw[,,R/keep],
                    out_oldtaus = taudraw[,,R/keep],
                    out_oldGamma = matrix(Gammadraw[R/keep,],byrow=TRUE,ncol=(nvars*npcov)),
                    out_oldVbeta = matrix(Vbetadraw[R/keep,],byrow=TRUE,ncol=nvars),
                    out_oldtheta = matrix(thetadraw[R/keep,],byrow=TRUE,ncol=nlvls),
                    out_bstep = bstep)
      }
      if (het_ind==1) {
        Cont = list(out_oldbetas = betadraw[,,R/keep],
                    out_oldtaus = taudraw[,,R/keep],
                    out_oldtheta = thetadraw[,,R/keep],
                    out_oldGamma = matrix(Gammadraw[R/keep,],byrow=TRUE,ncol=(nvars*npcov)),
                    out_oldVbeta = matrix(Vbetadraw[R/keep,],byrow=TRUE,ncol=nvars),
                    out_oldOmega = matrix(Omegadraw[R/keep,],byrow=TRUE,ncol=((nscov+1)*nlvls)),
                    out_bstep = bstep,out_ostep = ostep)
      }
    }
  }
  
  # Print total run time.
  ctime = proc.time()[3]
  cat(" Total Time Elapsed (in Hours/Minutes): ",round(((ctime - itime)/60)/60,2),"/",round((ctime - itime)/60,2),fill = TRUE)
  
  # Output.
  if (het_ind==0) {
    return(list(betadraw=betadraw,taudraw=taudraw,Gammadraw=Gammadraw,Vbetadraw=Vbetadraw,thetadraw=thetadraw,
                llikedraw=llikedraw,baccept=baccept,bstepdraw=bstepdraw,Cont=Cont))
  }
  if (het_ind==1) {
    return(list(betadraw=betadraw,taudraw=taudraw,thetadraw=thetadraw,Gammadraw=Gammadraw,Vbetadraw=Vbetadraw,Omegadraw=Omegadraw,
                llikedraw=llikedraw,baccept=baccept,oaccept=oaccept,bstepdraw=bstepdraw,ostepdraw=ostepdraw,Cont=Cont))
  }
}