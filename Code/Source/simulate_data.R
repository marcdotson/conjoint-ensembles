simulate_data <- function(
  nhh = 100,     # Number of respondents (households).
  ntask = 12,    # Number of choice tasks.
  nalt = 3,      # Number of choice alternatives.
  natt = 5,      # Number of attributes.
  nlevel = 3,    # Number of estimable attribute levels for each attribute.
  nversion = 10, # Number of versions of the experimental design.
  ana = FALSE,   # Attribute non-attendance flag.
  screen = FALSE # Screening pathology flag.
) {

  # Function to simulate choice data with or without pathologies.
  require(AlgDesign)
  
  # Create a full factorial design design and sample from it.
  data <- gen.factorial(rep(nlevel, times = natt), ntask, center = FALSE)
  indvec <- c(1:nrow(data))
  X.des <- matrix(double(ntask * nalt * nversion * (ncol(data) + 3)), ncol = (ncol(data) + 3))
  ii <- 1
  for(i in 1:nversion) {
    for(j in 1:ntask) {
      ind <- sample(indvec, nalt)
      xtemp <- data[ind,]
      xtemp <- cbind(c(1:nalt), xtemp)
      X.des[c(ii:(ii + (nalt - 1))), 1] <- i
      X.des[c(ii:(ii + (nalt - 1))), 2] <- j
      X.des[c(ii:(ii + (nalt - 1))), -c(1:2)] <- as.matrix(xtemp)
      ii <- ii + nalt
    }
  }
  
  # Create a dummy-coded design matrix.
  data.tmp <- as.data.frame(X.des)
  for(i in 4:ncol(data.tmp)) {
    data.tmp[,i] <- as.factor(data.tmp[,i])
  }
  data.tmp$y <- matrix(rnorm(nrow(data.tmp)), ncol = 1)
  out = lm(y ~ ., data.tmp, x = TRUE)
  desmat = out$x[,-1]
  
  # Generate betas and impose pathologies if flagged.
  nbeta <- natt * nlevel - natt
  ana.vec <- rep(1, times = nbeta)
  screen.vec <- rep(0, times = nbeta)
  if(ana == TRUE) {
    ana.draw <- sample(1:nbeta, 3, replace = FALSE)
    ana.vec[ana.draw] <- 0
  }
  if(screen == TRUE) {
    screen.draw <- sample(1:nbeta, 1, replace = FALSE)
    screen.vec[screen.draw] <- -100
  }
  bbar <- runif(nbeta, -1, 2)
  
  # Generate Y|X.
  regdata <- NULL
  for(i in 1:nhh) {
    X <- array(double(ntask * nalt * nbeta), dim = c(nalt, nbeta, ntask))
    Y <- double(ntask)
    nver <- sample(1:nversion, 1)
    tmp <- desmat[which(desmat[,1] == nver),]
    for(j in 1:ntask) {
      xtmp <- as.matrix(tmp[which(tmp[,2] == j), 4:(ncol(desmat))])
      betah <- bbar + rnorm(length(bbar), 0, 1)
      betah <- betah * ana.vec + screen.vec
      U <- as.vector(xtmp %*% betah) - as.vector(log(-log(runif(nalt))))
      Y[j] = which(U == max(U))
      X[,,j] = xtmp
    }
    regdata[[i]] = list(X = X, Y = Y)
  }
  
  # Format data for Stan.
  X <- array(double(nhh * ntask * nalt * nbeta), dim = c(nhh, ntask, nalt, nbeta))
  Y <- matrix(double(nhh * ntask), ncol = ntask)
  for(i in 1:nhh) {
    Y[i,] <- as.vector(regdata[[i]]$Y)
    for(j in 1:ntask) {
      X[i,j,,] <- regdata[[i]]$X[,,j]
    }
  }
  
  return(
    list(
      X = X,      # Design matrix.
      Y = Y,      # Choice data.
      bbar = bbar # Average betas.
    )
  )
}
