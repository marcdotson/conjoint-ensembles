clever_randomization <- function(
  Y = NULL,        # Choice data to cleverly randomize.
  X = NULL,        # Design matrices to cleverly randomize.
  natt = NULL,     # Number of attributes across design matrices.
  nlevels = NULL,  # Vector of number of attribute levels for each attribute.
  pct_test = .20,  # Percent of data to be saved for testing.
  nmember = 100,   # Number of members in the ensemble.
  test = ind_test, # Test flag.
  data = data      # Simulated data list.
) {
  # This function splits and returns training and testing data and 
  # two matrices that indicates if an attribute level is used for 
  # screening or is not attended to.
  
  # Needs to be modified to all for control over the type of pathology investigated 
  # and the number of attribute levels used for both screening and ANA.
  
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
      # Randomize respondent quality such that we continue with the same number of bootstrapped respondents.
      mat_resp[i,] <- sort(sample(1:nobs_train, nobs_train, replace = TRUE))
    }
    
    return(list(train_Y = train_Y, train_X = train_X, test_Y = test_Y,
                test_X = test_X, mat_ana = mat_ana, mat_screen = mat_screen,
                mat_resp = mat_resp))
  }
  if (test == 1) {
    mat_ana <- mat_screen <- matrix(double(nmember * nlevels_tot), nrow = nmember)
    nobs_train <- nrow(train_Y)
    mat_resp <- matrix(double(nmember * nobs_train), ncol = nobs_train)
    for (i in 1:nmember) {
      mat_ana[i,] <- data$mat_ana[1,]
      mat_screen[i,] <- data$mat_screen[1,]
      mat_resp[i,] <- data$mat_resp[train_ind]
    }
    
    return(list(train_Y = train_Y, train_X = train_X, test_Y = test_Y,
                test_X = test_X, mat_ana = mat_ana, mat_screen = mat_screen,
                mat_resp = mat_resp, betas = data$betas, bbar = data$bbar))
  }
}
