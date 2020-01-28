#********************************************************
#
#   Independence Metropolis algorithm for a HMNL model 
#                           with attribute non-attendance
#
#********************************************************
ana_hmnl=function(Data, Prior, Mcmc, TrueV){
  # This function implements a random-walk Metropolis-Hastings algorithm for a hierarchical
  # MNL with a multivariate normal distribution of heterogeneity and attribute non attendance
  
  # Data - list(nalt, data, lvlvec, Z)
  data=Data$data
  nalt=Data$nalt
  
  #MCMC values
  R = Mcmc$R  # Number of iterations in the Markov chain
  keep = Mcmc$keep  # Thinning parameter
  print = Mcmc$print  # Skip interval for printing output
  step = Mcmc$step  # Step size
  sim_ind=Mcmc$sim_ind
  out_ind=Mcmc$out_ind
  
  # Get argument values
  ntask=length(data[[1]]$y)/nalt  # Number of task in each scenario
  nresp=length(data)  # Total number of respondents
  lvlvec=Data$lvlvec-1 # Number of levels for each attribute (after removing null level)
  natt=length(lvlvec)  # Total number of attributes
  levloc=NULL 
  for(att in 1:natt){
    if(att==1){
    levloc[[att]]=list(loc=c(1:(lvlvec[att])))
    }else{ levloc[[att]]=list(loc=(cumsum(lvlvec)[att-1]+1):(cumsum(lvlvec)[att]))}
  }
  nbeta=ncol(data[[1]]$X)
  nz=ncol(Data$Z)
  Z=Data$Z
  
  # Get true values if they exist
  if(sim_ind==1){
    tDelta=TrueV$tDelta
    ttheta=TrueV$ttheta
    tbetamat=t(TrueV$tbetamat)
    tC=t(TrueV$tC)
  }

  # Prior = list(Deltabar,ADelta, nu, V, alpha).
  Deltabar = Prior$Deltabar  #Means for upper level parameters
  ADelta = Prior$ADelta  #Precision matrix for normal prior on Delta.
  nu = Prior$nu  #DF for IW prior on Vbeta.
  V = Prior$V  #Location for IW prior on Vbeta.
  alpha = Prior$alpha  #hyperparameters for beta prior on tau
  
  #Starting values
  oldbetamat=matrix(double(nbeta*nresp),nr=nresp)  #initial betas (partworths)
  oldbetastarmat=matrix(double(nbeta*nresp),nr=nresp)
  oldCmat= matrix(double(nbeta*nresp)+1,nr=nresp)  #initial taus (attendance values)
  oldtheta=c(rep(.5,natt))  #initial prob for attendance
  oldDelta=matrix(double(nz*nbeta),nr=nz)  #initial value for Upper level coefficients
  oldbetabarmat=Z%*%oldDelta  #initial betabar (mean)
  oldVbeta=.5+diag(nbeta)*.5  #initial Vbeta
  oldVbetai=backsolve(chol(oldVbeta),diag(nbeta))  #initial inv chol Vbeta
  llike=matrix(double(nresp),nc=1) #log likelihood
  
  #Setup functions
  #function that returns the value of the prior for betas
  logprior=function(beta,betabar,Vbi){
    return((beta-betabar)%*%Vbi%*%t(beta-betabar)*(-.5))
  }
  
  McondMom=function(x,mu,sig,i){
    sigii = sig[i,i]
    sigini = sig[i,-i]
    signini = sig[-i,-i]
    signii = sig[-i,i]
    m = mu[i] +  sigini%*%chol2inv(chol(signini))%*%(x[-i] - mu[-i])
    s = sigii - sigini%*%chol2inv(chol(signini))%*%signii
    return(list(cmean = as.vector(m), cvar = s))
  }
    
  #loglike function
  loglike=function(xb,y){
    probs=log(exp(xb))-log(matrix(rep(colSums(exp(xb)),length(
      xb[,1])),byrow=T,nc=length(xb[1,])))
    loc=cbind(y,c(1:(ncol(xb))))
    return(sum(probs[loc]))
  }

  #Setup storage
  betadraws=array(double(nbeta*nresp*R/keep),dim=c(R/keep,nresp,nbeta))
  betastardraws=array(double(nbeta*nresp*R/keep),dim=c(R/keep,nresp,nbeta))
  Deltadraws=array(double(nbeta*nz*R/keep),dim=c(R/keep,nbeta,nz))
  Vbetadraws=array(double(nbeta*nbeta*R/keep),dim=c(R/keep,nbeta,nbeta))
  Cdraws=array(double(nbeta*nz*R/keep),dim=c(R/keep,nresp,nbeta))
  thetadraws=matrix(double(natt*R/keep),nc=natt)
  llikedraw=matrix(double(nresp*R/keep),nr=R/keep)

  #timer
  itime = proc.time()[3]
  
  for(r in 1:R){
    oldbetabarmat=Z%*%oldDelta
    attcount=matrix(double(nresp*natt),nc=natt)
    #Draw new betas and attendance parameters for each respondent
    for(resp in 1:nresp){
      oldbeta=oldbetamat[resp,]
      oldbetabar=oldbetabarmat[resp,]
      oldC=oldCmat[resp,]
      for(j in 1:natt){
        #draw each tau and beta together
        newtj=rbinom(1,1,oldtheta[j]) #draw tau for this attribute/respondent
        newtj=newtj-(newtj-1)*.001 #cannot be zero
        i=levloc[[j]]$loc  #identify locations of betas for this attribute
        newC=oldC;newC[i]=newtj  #create new C
        moms=McondMom(oldbeta,diag(newC)%*%oldbetabar,
                      diag(newC)%*%oldVbeta%*%diag(newC),i) #get conditional moments
        newbeta=oldbeta  #setup new beta draw
        newbeta[i]=t(chol(moms$cvar))%*%rnorm(lvlvec[j])+moms$cmean  #draw from cond mult norm

        #get likelihood and prior values
        oldllikeb=loglike(matrix(data[[resp]]$X%*%oldbeta,nr=nalt)
                          ,as.matrix(data[[resp]]$y))
        newllikeb=loglike(matrix(data[[resp]]$X%*%newbeta,nr=nalt)
                          ,as.matrix(data[[resp]]$y))
        diffvecb=newllikeb-oldllikeb
        alphab=min(exp(diffvecb), 1)
      
        #accept or reject
        draw=runif(1)
        acceptb=0
        if(alphab>draw){acceptb=1}
          llike[resp]=oldllikeb
          if(r>10){
        if(acceptb==1){
          oldC=newC
          oldbeta=newbeta
          llike[resp]=newllikeb
        }}
        attcount[resp,j]=round(oldC[levloc[[j]]$loc][1],1) #store attendance for this resp/att
      }
      oldbetamat[resp,]=oldbeta
      oldbetastarmat[resp,]=chol2inv(chol(diag(oldC)))%*%oldbeta
      oldCmat[resp,]=oldC
    }
  
    #Draw new values for Delta
    outbetaup=rmultireg(oldbetastarmat,matrix(Z,nc=nz),matrix(Deltabar,nr=nz),ADelta,nu,V)
    oldDelta=outbetaup$B
    oldVbeta=outbetaup$Sigma
    oldVbetai=chol2inv(chol(oldVbeta))
    
    #Draw new values for theta
    for(j in 1:(natt)){
      oldtheta[j]=rbeta(1,alpha[1]+sum(attcount[,j]),alpha[1]+sum(1-attcount[,j]))
    }
    
    #Store values
    if(r%%keep==0){
      betadraws[r/keep,,]=oldbetamat
      betastardraws[r/keep,,]=oldbetastarmat
      Deltadraws[r/keep,,]=oldDelta
      Vbetadraws[r/keep,,]=oldVbeta
      Cdraws[r/keep,,]=oldCmat
      thetadraws[r/keep,]=oldtheta
      llikedraw[r/keep,]=llike
    }
    
    #print progress
    # Print progress.
    if (r%%(keep*5) == 0){
      ctime = proc.time()[3]
      timetoend = ((ctime - itime)/r)*(R - r)
      cat(" ",r," ( hrs/min to end",round((timetoend/60)/60,2),"/",round(timetoend/60,2),"| step",
            round(step,5),"| ll",round(sum(llike),2),")",fill = TRUE)
    }

    
    #print chart
    if(r%%print==0){
      par(mfrow=c(3,1))
      plot(rowSums(llikedraw),type="l",ylab="Log Likelihood")
      matplot(matrix(Deltadraws,nr=R/keep),type="l",ylab="Deltabar Draws")
      if(sim_ind==1){abline(h=tDelta)}
      matplot(thetadraws,type="l",ylab="theta Draws")
      if(sim_ind==1){abline(h=ttheta)}
    }
  } 
  return(list(betadraws=betadraws,betastardraws=betastardraws,Deltadraws=Deltadraws,  
              Cdraws=Cdraws,thetadraws=thetadraws,Vbetadraws=Vbetadraws,llikes=llikes))
}
 
  
 