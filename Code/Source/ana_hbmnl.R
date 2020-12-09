## Esimation routine for attribute non-attendence using the Gilbride, Allenby and Brazell (2006)
## JMR model for heterogeneous variable selection.  Variable selection occurs at the attribute-level
## Required data structure for oberservations follows the conventions of bayesm.  Data is a list of lists contain
## X (design matrix) & Y (multinomial choice indicator) objects

rana_hmnl <- function(Data,Z,R=1000,keep=1,ss.b=.4, amap = NULL){
  
  ## load libraries
  library(bayesm)
  
  ################################################################################
  
  ## internal functions
  
  rmultiregG = function (Y, X, Bbar, A, nu, V) {
    n = nrow(Y)
    m = ncol(Y)
    k = ncol(X)
    RA = chol(A)
    W = rbind(X, RA)
    Z = rbind(Y, RA %*% Bbar)
    IR = backsolve(chol(crossprod(W)), diag(k))
    Btilde = crossprod(t(IR)) %*% crossprod(W, Z)
    S = crossprod(Z - W %*% Btilde)
    rwout = rwishart(nu + n, chol2inv(chol(V + S)))
    B = Btilde + IR %*% matrix(rnorm(m * k), ncol = m) %*% t(rwout$CI)
    return(list(B = B, Sigma = rwout$IW))
  }
  ################################################################################
  
  ## Unpack Data
  X = Data[[1]]$X[,,1] # Exemplar design matrix: array(K,A,J)
  Y = Data[[1]]$Y  # J*1 matrix of multinomial choices
  
  ################################################################################
  
  ## define constants
  J = length(Y) #number of choice taskts
  K1 = nrow(X) #number of alternatives per choice set
  A = ncol(X) #total number of attribute levels
  N = length(Data) #number of respondents
  nZ = ncol(Z)
  Sig1 = IMat1 = diag(K1)
  sigi1=chol2inv(chol(Sig1))
  
  ss.b = .07 ## initial step size for m-h step.  Will be modified dynamically at the individual-level
  ssB = double(N) + ss.b
  
  ## amat defines the mapping between attribute levels in the design matrix and their corresponding attribute indicators
  ## the default selection will lead to attribute-level non-attendence. 
  if(is.null(amat)) amat = matrix(c(1:A,1:A),ncol=2)
  nM = max(amat[,2]) ## Number of attribtuts as defined by the amat mapping
  
  ################################################################################
  
  ## allocate storage space
  betadraw = array(double((R/keep)*N*A),dim=c(R/keep,A,N))
  bbardraw = matrix(double((R/keep)*A),ncol=A)
  vbetakeep = matrix(double(A*A*R/keep),ncol=A*A)
  nacceptB = matrix(double(N),ncol=1)
  lltemp = matrix(double(N),ncol=1)
  llkeep = matrix(double(N*R/keep),ncol=N)
  thetakeep = matrix(double(nM,(R/keep)),ncol=nM)

  ################################################################################
  ## Starting values for beta & theta
  
  ## START HERE ##
  
  ## Uniform [-.01,.01]
  betaM = betaMnew = matrix(runif(N*A,-.01,.01),ncol=A) 
  theta = double(nM) + .5
  
  ## Initialize variables
  beta = c(double(A))
  bbar = c(double(A))
  ztemp = c(1:K1)
  LB = -50
  UB = 50
  iota = matrix(1,ncol=1,nrow=N)
  
  ## Prior values 
  bbar.p = c(double(A)) + 2
  nu = A+3
  A.b = 0.01
  V.b = Sig.b = diag(A)*nu #*.1 pp 74 in text
  rooti.b = backsolve(chol(V.b*5),diag(A))
  
  ## mcmc variables starting values
  llold.b = 0
  pold.b = lndMvn(beta,bbar,rooti.b)
  
  ################################################################################
  
  ## Estimation Algorithm
  
  cat(" ", fill = TRUE)
  cat("Starting Estimation Routine for Hierarchical MNL:", fill = TRUE)
  itime = proc.time()[3]
  cat("MCMC Iteration (est time to end - min) ", fill = TRUE)
  fsh()
  for(rep in 1:R){
    
    ## Unit Loop 
    for(nhh in 1:N){
      
      ################################################################################
      
      ## unpack data for household nhh
      X1.nhh = Data[[nhh]]$X
      Y1 = Data[[nhh]]$Y
      betaold = betaM[nhh,]
      
      llold = lltemp[nhh,]
      ss.b = ssB[nhh]
      
      ################################################################################
      ## Step 1. Draw beta|else - via M-H (evaluate likelihood using GHK)
      
      if(rep==1){
        ## Compute current value of likelihood & prior
        llold = 0
        
        ## Compute for first set of choices (4 alt)
        for(ntask in 1:J){
          Y = Data[[nhh]]$Y[ntask]
          X = X1.nhh[,,ntask]
          mu = X%*%betaold
          ltemp = exp(mu)/sum(exp(mu))
          llold = llold + log(ltemp[Y])
        }
        
      }
      
      pold.b = lndMvn(betaold,as.vector(bbar),rooti.b)
      
      ## draw betanew
      betanew = betaold + rnorm(length(betaold),0,ss.b)
      
      ## Compute likelihood & prior for new values of beta
      llnew = 0  
      for(ntask in 1:J){
        Y = Data[[nhh]]$Y[ntask]
        X = X1.nhh[,,ntask]
        mu = X%*%betanew
        ltemp = exp(mu)/sum(exp(mu))
        llnew = llnew + log(ltemp[Y])
      }
      
      pnew.b = lndMvn(betanew,as.vector(bbar),rooti.b)
      
      ## metropolis step
      ldiff = llnew + pnew.b - llold - pold.b
      
      alpha = min(1,exp(ldiff))
      if (alpha<1) {unif=runif(1)} else {unif=0}
      
      ## update variables and save draws if accepted
      if (unif <= alpha)
      {
        llold = llnew
        betaold = betanew
        nacceptB[nhh,] = nacceptB[nhh,] + 1 } else { }
      
      betaM[nhh,] = betaold
      betaMnew[nhh,] = betanew
      
      lltemp[nhh,] = llold
      
    } ## end unit loop
    
    ################################################################################
    ## Step 4. Draw hierarchical priors
    if(N>1){
      
      ## Distribution of heterogeneity for beta
      rmout = rmultiregG(betaM,iota,bbar.p,A.b,nu,V.b)
      bbar = rmout$B
      Sig.b = rmout$Sigma
      rooti.b = backsolve(chol(rmout$Sigma),diag(A))
      
    } ## End N>1 condition
    
    ################################################################################
    ## Print Time & Parameter Estimates
    if (rep%%10 == 0) {
      ctime = proc.time()[3]
      timetoend = ((ctime-itime)/rep) * (R-rep)
      cat(" ", rep, " ", round(timetoend/60,1), " ", round(bbar,2),
          " ", round(mean(nacceptB/10),2)," ",
          sum(lltemp), " ", fill = TRUE)
      fsh()
    }
    
    ################################################################################
    ## Save Draws
    if (rep%%keep == 0){
      betadraw[rep/keep,,] = t(betaM) 
      bbardraw[rep/keep,] = bbar
      vbetakeep[rep/keep,] = as.vector(rmout$Sigma)
      llkeep[rep/keep,] = as.vector(lltemp)
    }
    
    ## Dynamically update step size
    
    if(rep%%10 == 0){
      if(rep<R/2){
        for(ii in 1:N){
          if(nacceptB[ii]/10 < .2) ssB[ii] = ssB[ii]*.9
          if(nacceptB[ii]/10 > .4) ssB[ii] = ssB[ii]*1.1
        }
      }
      nacceptB = nacceptB*0
    }
    
  } ## end MCMC
  
  return(list(betadraw=betadraw, bbardraw=bbardraw, vbetadraw=vbetakeep,
              llkeep = llkeep))
}


