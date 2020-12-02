sim_data <- function(nhh = 100, nalt = 3, ntask = 12, natt = 5, nlevel = 3, 
                     nversion = 10, ana = FALSE, screen = FALSE){

  ## Function to simulate data with or without pathologies
  ## nhh = nubmer of households; nalt = number of alternatives per choice; ntask = number of tasks
  ## natt = number attributes; nlevel = number of levels per attribute (consider making variable)
  ## ana = attribute non-attendance flag; screen = screening flag
  require(AlgDesign)
  
  # nhh = 100
  # nalt = 3
  # ntask = 12 
  # natt = 5 
  # nversion = 10
  # ana = FALSE 
  # screen = FALSE
  # nlevel = 3 
    
  dat <- gen.factorial(rep(nlevel,times=natt),ntask,center=FALSE)

  ## create a design by sampling from the full factorial design
  
  indvec = c(1:nrow(dat))
  
  X.des = matrix(double(ntask*nalt*nversion*(ncol(dat)+3)),ncol=(ncol(dat)+3))
  
  ii = 1
  for(i in 1:nversion){
    for(j in 1:ntask){
      ind = sample(indvec,nalt)
      xtemp = dat[ind,]
      xtemp = cbind(c(1:nalt),xtemp)
      X.des[c(ii:(ii+(nalt-1))),1] = i
      X.des[c(ii:(ii+(nalt-1))),2] = j
      X.des[c(ii:(ii+(nalt-1))),-c(1:2)] = as.matrix(xtemp)
      ii = ii + nalt
    }
  }
  
  ## Create Dummy Coded Design Matrix
  dat.tmp = as.data.frame(X.des)
  
  for(i in 4:ncol(dat.tmp)){
    dat.tmp[,i] = as.factor(dat.tmp[,i])
  }
  
  dat.tmp$y = matrix(rnorm(nrow(dat.tmp)),ncol=1)
  out = lm(y~.,dat.tmp,x=TRUE)
  desmat = out$x[,-1]
  
  ## Generate betas, Y|X, format data for Stan
  nbeta = natt*nlevel - natt
  
  if(ana==TRUE){
    ana.vec = rep(1,times = nbeta)
    ana.draw = sample(1:nbeta, 3, replace = FALSE)
    ana.vec[ana.draw] = 0
  }

  if(screen==TRUE){
    screen.vec = rep(0,times = nbeta)
    screen.draw = sample(1:nbeta, 1, replace = FALSE)
    screen.vec[screen.draw] = -100
  }
  
    
  bbar = runif(nbeta,-1,2)
  
  regdata = NULL
  
  for(i in 1:nhh){
    X = array(double(ntask*nalt*nbeta),dim=c(nalt,nbeta,ntask))
    Y = double(ntask)
    nver = sample(1:nversion,1)
    tmp = desmat[which(desmat[,1]==nver),]
    for(j in 1:ntask){
      xtmp = as.matrix(tmp[which(tmp[,2]==j),4:(ncol(desmat))])
      betah = bbar + rnorm(length(bbar),0,1)
      betah = betah * ana.vec + screen.vec
      U = as.vector(xtmp%*%betah) - as.vector(log(-log(runif(nalt))))
      Y[j] = which(U==max(U))
      X[,,j] = xtmp
    }
    regdata[[i]] = list(X = X, Y = Y)
  }
  
  ##Format data for Stan
  X = array(double(nhh*ntask*nalt*nbeta),dim=c(nhh,ntask,nalt,nbeta))
  Y = matrix(double(nhh*ntask),ncol=ntask)
  
  for(i in 1:nhh){
    Y[i,] = as.vector(regdata[[i]]$Y)
    for(j in 1:ntask){
      X[i,j,,] = regdata[[i]]$X[,,j]
    }
  }
  
  return(list(X=X, Y=Y, bbar = bbar))
  
}