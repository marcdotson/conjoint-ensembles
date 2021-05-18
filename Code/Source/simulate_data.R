simulate_data <- function(
  nhh = 100,      # Number of households/respondents.
  ntask = 12,     # Number of choice tasks.
  nalt = 3,       # Number of choice alternatives.
  natt = 5,       # Number of discrete attributes.
  nlevel = 3,     # Number of levels for each discrete attribute.
  nversion = 10,  # Number of versions of the experimental design.
  ana = FALSE,    # Attribute non-attendance flag.
  screen = FALSE, # Screening pathology flag.
  resp = FALSE,   # Respondent quality flag.
  hetero = FALSE  # Pathologies at the individual-level.
) {

  # Function to simulate choice data with or without pathologies.
  require(AlgDesign)
  
  # Create a full factorial design design and sample from it.
  data <- gen.factorial(rep(nlevel, times = natt), ntask, center = FALSE)
  indvec <- c(1:nrow(data))
  X.des <- matrix(double(ntask * nalt * nversion * (ncol(data) + 3)), ncol = (ncol(data) + 3))
  ii <- 1
  for (i in 1:nversion) {
    for (j in 1:ntask) {
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
  for (i in 4:ncol(data.tmp)) {
    data.tmp[,i] <- as.factor(data.tmp[,i])
  }
  data.tmp$y <- matrix(rnorm(nrow(data.tmp)), ncol = 1)
  out = lm(y ~ ., data.tmp, x = TRUE)
  desmat = out$x[,-1]
  
  # Generate betas and impose pathologies if flagged.
  nbeta <- natt * nlevel - natt
  # ana.vec <- rep(1, times = nbeta)
  # screen.vec <- rep(0, times = nbeta)
  ana.mat <- matrix(double(nbeta * nhh), ncol = nbeta) + 1
  screen.mat <- matrix(double(nbeta * nhh), ncol = nbeta)
  
  for (i in 1:nhh) {
    if (ana == TRUE) {
      # ana.draw <- sample(1:nbeta, 3, replace = FALSE)
      # ana.mat[i,ana.draw] <- 0
      ana.draw <- sample(1:natt, size = round(runif(n = 1, min = 1, max = natt - 1)), replace = FALSE)
      for (j in 1:natt) {
        if (j %in% ana.draw) ana.mat[i, ((j * nlevel - j) - 1):(j * nlevel - j)] <- 0
      }
    }
    if (screen == TRUE) {
      screen.draw <- sample(1:nbeta, 1, replace = FALSE)
      screen.mat[i, screen.draw] <- -100
    }
  }
  
  # Induce randomization in simulated data for respondent data 
  # forcing random choice if resp == TRUE.
  resp.id <- 0
  if (resp == TRUE) {
    resp.id <- sample(1:nhh, sample(0:25, 1), replace = FALSE)
  }
  
  if (hetero == FALSE) {
    ana.vec <- ana.mat[1,]
    screen.vec <- screen.mat[1,]
    for (i in 1:nhh) {
      ana.mat[i,] <- ana.vec
      screen.mat[i,] <- screen.vec
    }
  }
  
  bbar <- runif(nbeta, -1, 2)
  
  # Generate Y|X.
  regdata <- NULL
  for (i in 1:nhh) {
    X <- array(double(ntask * nalt * nbeta), dim = c(nalt, nbeta, ntask))
    Y <- double(ntask)
    nver <- sample(1:nversion, 1)
    tmp <- desmat[which(desmat[,1] == nver),]
    ana.vec <- ana.mat[i,]
    screen.vec <- screen.mat[i,]
    for (j in 1:ntask) {
      xtmp <- as.matrix(tmp[which(tmp[,2] == j), 4:(ncol(desmat))])
      betah <- bbar + rnorm(length(bbar), 0, 1)
      betah <- betah * ana.vec + screen.vec
      U <- as.vector(xtmp %*% betah)
      prob.y = exp(U)/sum(exp(U))
      if(i %in% resp.id) prob.y = prob.y*0 + 1/nalt
      Y[j] = sample(1:nalt,1,prob = prob.y)
      #U <- as.vector(xtmp %*% betah) - as.vector(log(-log(runif(nalt)))) + resp.err
      #Y[j] = which(U == max(U))
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
