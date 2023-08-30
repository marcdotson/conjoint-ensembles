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
  hetero = FALSE, # Pathologies at the individual-level.
  test = FALSE    # Test flag.
) {

  # Function to simulate choice data with or without pathologies.
  require(AlgDesign)
  
  # Create a full factorial design and sample from it where the design only has 
  # discrete attributes and each attribute has the same number of levels.
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
  
  # Generate betas and prepare to impose pathologies if flagged.
  nbeta <- natt * nlevel - natt
  ana.mat <- matrix(double(nbeta * nhh), ncol = nbeta) + 1
  screen.mat <- matrix(double(nbeta * nhh), ncol = nbeta)
  
  for (i in 1:nhh) {
    # If ANA is flagged, simulate ANA where each respondent pays attention to at least one attribute and 
    # has non-attendance for at least one attribute.
    if (ana == TRUE) {
      ana.draw <- sample(1:natt, size = round(runif(n = 1, min = 1, max = natt - 1)), replace = FALSE)
      for (j in 1:natt) {
        if (j %in% ana.draw) ana.mat[i, ((j * nlevel - j) - 1):(j * nlevel - j)] <- 0
      }
    }
    # If screening is flagged, simulate screening where each respondent screens based on at least one 
    # attribute level but not all of them.
    if (screen == TRUE) {
      screen.draw <- sample(1:nbeta, size = round(runif(n = 1, min = 1, max = nbeta -1)), replace = FALSE)
      screen.mat[i, screen.draw] <- -100
    }
  }
  
  # If respondent quality is flagged, randomly determine up to 25 low-quality respondents.
  # force a random choice for up to 25 respondents.
  ####################################
  # REWRITE TO MIRROR ANA & SCREENING (always heterogeneous).
  resp.id <- 0
  if (resp == TRUE) {
    resp.id <- sample(1:nhh, sample(0:25, 1), replace = FALSE)
  }
  ####################################
  
  # If heterogeneity isn't flagged, have the first iteration of ANA and screening apply for all members.
  if (hetero == FALSE) {
    ana.vec <- ana.mat[1,]
    screen.vec <- screen.mat[1,]
    for (i in 1:nhh) {
      ana.mat[i,] <- ana.vec
      screen.mat[i,] <- screen.vec
    }
  }
  
  # Generate average betas.
  bbar <- runif(nbeta, -1, 2)
  betas <- NULL
  
  # Generate Y|X.
  regdata <- NULL
  for (i in 1:nhh) {
    X <- array(double(ntask * nalt * nbeta), dim = c(nalt, nbeta, ntask))
    Y <- double(ntask)
    nver <- sample(1:nversion, 1)
    tmp <- desmat[which(desmat[,1] == nver),]
    ana.vec <- ana.mat[i,]
    screen.vec <- screen.mat[i,]
    # Generate betas as a deviation from the average.
    betah <- bbar + rnorm(length(bbar), 0, 1)
    betas[[i]] <- betah
    for (j in 1:ntask) {
      # Multiply betas by an ANA indicator and add screening (defaults to 1 and 0, respectively).
      xtmp <- as.matrix(tmp[which(tmp[,2] == j), 4:(ncol(desmat))])
      betah <- betas[[i]] * ana.vec + screen.vec
      U <- as.vector(xtmp %*% betah)
      prob.y <- exp(U) / sum(exp(U))
      # If respondent quality is flagged, force a random choice for the up to 25 respondents.
      if (i %in% resp.id) prob.y <- prob.y * 0 + 1 / nalt
      Y[j] <- sample(1:nalt, 1, prob = prob.y)
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
  
  if (test == 0) {
    return(
      list(
        X = X,      # Design matrix.
        Y = Y,      # Choice data.
        bbar = bbar # Average betas.
      )
    )
  }
  if (test == 1) {
    mat_ana <- ifelse(ana.mat == 0, 1, 0)
    mat_screen <- ifelse(screen.mat != 0, 1, 0)
    ####################################
    # REWRITE TO MIRROR ANA & SCREENING (always heterogeneous).
    # mat_resp <- 1:nhh
    # mat_resp <- sort(c(mat_resp[-resp.id], mat_resp)[1:nhh])
    ####################################
    
    return(
      list(
        X = X,      # Design matrix.
        Y = Y,      # Choice data.
        mat_ana = mat_ana,
        mat_screen = mat_screen,
        # mat_resp = mat_resp,
        bbar = bbar, # Average betas.
        betas = betas
      )
    )
  }
  
}
