# Preamble ----------------------------------------------------------------
# Load packages.
library(tidyverse)

# Set the simulation seed.
set.seed(42)

# Simulate Data and Induce Clever Randomization ---------------------------
if (ind_sim == 1) {
  if (!file.exists(here::here("data", str_c("sim_", file_id, ".rds")))) {
    # Simulate data, induce clever randomization, and save.
    data <- simulate_data(
      nhh = 300,           # Number of respondents (households).
      ntask = 12,          # Number of choice tasks.
      nalt = 3,            # Number of choice alternatives.
      natt = 5,            # Number of attributes.
      nlevel = 3,          # Number of attribute levels for each attribute.
      nversion = 10,       # Number of versions of the experimental design.
      ana = ind_ana,       # Attribute non-attendance flag.
      screen = ind_screen, # Screening flag.
      resp = ind_resp,     # Respondent quality flag.
      hetero = ind_hetero, # Pathologies differ by individual flag.
      test = ind_test      # Test flag.
    )
    
    ####################################################
    # Incorporate simulate_data.R
    ####################################################
    
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
    
    ####################################################
    # We don't include the probability of there being ANA if ana == TRUE.
    # We don't include the probability of there being screening if screen == TRUE.
    # Allows for screening on ANA levels as well as price.
    ####################################################
    
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
    # REWRITE TO MIRROR ANA & SCREENING (always heterogeneous) as a deviation from latent utility.
    resp.id <- 0
    if (resp == TRUE) {
      resp.id <- sample(1:nhh, sample(0:25, 1), replace = FALSE)
    }
    ####################################
    
    # If heterogeneity isn't flagged, have the first iteration of ANA and screening apply for all members.
    
    #######
    # Need to make this an array regardless of hetero.
    #######
    
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
          # mat_resp = mat_resp, # Requires REWRITE TO MIRROR ANA & SCREENING (always heterogeneous).
          bbar = bbar, # Average betas.
          betas = betas
        )
      )
    }
    
    ####################################################
    
    data <- clever_randomization(
      Y = data$Y,          # Choice data to cleverly randomize.
      X = data$X,          # Design matrices to cleverly randomize.
      natt = 5,            # Number of attributes across design matrices.
      nlevels = rep(3, 5), # Vector of number of attribute levels for each attribute.
      pct_test = .20,      # Percent of data to be saved for testing.
      nmember = 2000,      # Number of possible members in the ensemble.
      test = ind_test,     # Test flag.
      data = data          # Simulated data list.
    )
    
    ####################################################
    # Incorporate clever_randomization.R
    ####################################################
    
    Y = NULL,        # Choice data to cleverly randomize.
    X = NULL,        # Design matrices to cleverly randomize.
    natt = NULL,     # Number of attributes across design matrices.
    nlevels = NULL,  # Vector of number of attribute levels for each attribute.
    pct_test = .20,  # Percent of data to be saved for testing.
    nmember = 100,   # Number of members in the ensemble.
    test = ind_test, # Test flag.
    data = data      # Simulated data list.
    
    # This function splits and returns training and testing data and 
    # two matrices that indicates if an attribute level is used for 
    # screening or is not attended to.
    
    ####################################
    # Needs to be modified to all for control over the type of pathology investigated 
    # and the number of attribute levels used for both screening and ANA, for both
    # homogeneous and heterogeneous outcomes.
    ####################################
    
    ####################################
    # How do we think about covering the potential models space?
    ####################################
    
    # Determine the number of observation nunits and total number of levels.
    nobs <- dim(X)[1]
    nlevels_tot <- sum(nlevels) - natt
    
    # Randomly draw respondents in the testing data.
    ntest <- round(pct_test * nobs, 0)
    test_ind <- sort(sample(1:nobs, ntest, replace = FALSE))
    '%!in%' <- function(x, y)!('%in%'(x, y))
    train_ind <- which(c(1:nobs) %!in% test_ind)
    
    # Split into training and testing.
    train_Y <- Y[train_ind,]
    train_X <- X[train_ind,,,]
    test_Y <- Y[test_ind,]
    test_X <- X[test_ind,,,]
    
    if (test == 0) {
      # Create randomization patterns for each pathology, no flag necessary.
      mat_ana <- mat_screen <- matrix(double(nmember * nlevels_tot), nrow = nmember)
      nobs_train <- nrow(train_Y)
      mat_resp <- matrix(double(nmember * nobs_train), ncol = nobs_train)
      
      ####################################################
      # We don't include the probability of there being ANA if ana == TRUE.
      # We don't include the probability of there being screening if screen == TRUE.
      # Allows for screening on ANA levels as well as price.
      ####################################################
      
      for (i in 1:nmember) {
        # Randomize ANA such that each respondent pays attention to at least one attribute 
        # and has non-attendance for at least one attribute.
        ana.ind <- sample(1:natt, size = round(runif(n = 1, min = 1, max = natt - 1)), replace = FALSE)
        for (j in 1:natt) {
          if (j %in% ana.ind) mat_ana[i, ((sum(nlevels[1:j]) - j) - (nlevels[j] - 2)):(sum(nlevels[1:j]) - j)] <- 1
        }
        # Randomize screening such that each respondent screens based on at least one attribute level.
        screen.ind <- sample(1:nlevels_tot, size = round(runif(n = 1, min = 1, max = nlevels_tot)), replace = FALSE)
        mat_screen[i, screen.ind] <- 1
        # # Randomize respondent quality such that we continue with the same number of bootstrapped respondents.
        # mat_resp[i,] <- sort(sample(1:nobs_train, nobs_train, replace = TRUE))
      }
      
      return(list(train_Y = train_Y, train_X = train_X, test_Y = test_Y, test_X = test_X, 
                  mat_ana = mat_ana, mat_screen = mat_screen, mat_resp = mat_resp))
    }
    if (test == 1) {
      # Split pathology matrices into training.
      mat_ana <- data$mat_ana[train_ind,]
      mat_screen <- data$mat_screen[train_ind,]
      nobs_train <- nrow(train_Y)
      mat_resp <- matrix(double(nmember * nobs_train), ncol = nobs_train)
      
      return(list(train_Y = train_Y, train_X = train_X, test_Y = test_Y, test_X = test_X,
                  mat_ana = mat_ana, mat_screen = mat_screen, mat_resp = mat_resp,
                  betas = data$betas, bbar = data$bbar))
    }
    
    ####################################################
    
    # Save simulated data.
    write_rds(data, here::here("data", str_c("sim_", file_id, ".rds")))
  }
}

# Clean Data and Induce Clever Randomization ------------------------------
if (ind_emp == 1) {
  if (!file.exists(here::here("data", str_c("emp_", data_id, ".rds")))) {
    # Load design.
    design <- read.csv(here::here("data", str_c("emp_", data_id, "_design.csv")), header = TRUE)
    design <- design[,-1]
    
    # Load choice data.
    choice_data <- read.csv(here::here("data", str_c("emp_", data_id, "_final.csv")), header=TRUE)
    choice_data <- choice_data[,-1]
    
    # Sample Characteristics
    nobs <- nrow(choice_data)
    ntask <- 10
    nalt <- 4
    natt <- ncol(design) - 3
    
    # Format Data for Stan
    regdata <- NULL
    X <- array(double(nobs*ntask*nalt*natt),dim=c(nobs,ntask,nalt,natt))
    Y <- as.matrix(choice_data[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10")])
    
    for (i in 1:nobs) {
      nver <- choice_data[i,"X0.Version"]
      tmp <- design[which(design[,1]==nver),]
      for (j in 1:ntask) {
        X[i,j,,] <- as.matrix(tmp[which(tmp[,2] == j), 4:(ncol(design))])
      }
    }
    
    # # Induce clever randomization.
    # data <- clever_randomization(
    #   Y = Y,                                  # Choice data to cleverly randomize.
    #   X = X,                                  # Design matrices to cleverly randomize.
    #   natt = 9,                               # Number of attributes across design matrices.
    #   nlevels = c(5, 8, 3, 4, 2, 2, 2, 2, 2), # Vector of attribute levels for each attribute.
    #   pct_test = .1,                          # Percent of data to be saved for testing.
    #   nmember = 2000                          # Number of possible members in the ensemble.
    # )
    
    # Save data.
    write_rds(data, here::here("data", str_c("emp_", data_id, ".rds")))
  }
}

