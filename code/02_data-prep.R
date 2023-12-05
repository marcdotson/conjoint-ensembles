# Simulate Data and Induce Clever Randomization ---------------------------
if (ind_sim == 1) {
  if (!file.exists(here::here("data", str_c("sim_", file_id, ".rds")))) {
    # Generate a full factorial design and corresponding row index.
    design_full <- gen.factorial(rep(nlevel, times = natt), ntask, center = FALSE)
    index_full <- c(1:nrow(design_full))
    
    # Generate a fractional factorial design.
    design_fractional <- matrix(
      double(ntask * nalt * nversion * (ncol(design_full) + 3)), 
      ncol = (ncol(design_full) + 3)
    )
    first_alt <- 1
    for (version in 1:nversion) {
      for (task in 1:ntask) {
        # Sampling alternatives from the full factorial design and include version and
        # task indicators as part of each task in the fractional factorial design.
        index_task <- sample(index_full, nalt)
        design_task <- as.matrix(cbind(version, task, c(1:nalt), design_full[index_task,]))
        design_fractional[c(first_alt:(first_alt + (nalt - 1))), ] <- design_task
        
        # Increment the index for the first alt for the next task in the fractional design.
        first_alt <- first_alt + nalt
      }
    }
    
    # Generate a dummy-coded design matrix by (ab)using lm().
    design_dummy <- as.data.frame(design_fractional)
    for (level in 4:ncol(design_dummy)) {
      design_dummy[,level] <- as.factor(design_dummy[,level])
    }
    design_dummy$y <- matrix(rnorm(nrow(design_dummy)), ncol = 1)
    out <- lm(y ~ ., design_dummy, x = TRUE)
    design_dummy <- out$x[,-1]
    
    # Prepare to impose pathologies probabilistically, conditioned on indicator flags.
    nbeta <- natt * nlevel - natt
    mat_ana <- matrix(double(nbeta * nresp), ncol = nbeta) + 1
    mat_screen <- matrix(double(nbeta * nresp), ncol = nbeta)
    mat_qual <- matrix(double(nresp), ncol = 1)
    for (resp in 1:nresp) {
      # If ANA is flagged, with prob_ana, simulate ANA where a respondent pays attention to at 
      # least one attribute and has non-attendance for at least one attribute.
      if (ind_ana == TRUE & runif(1) < prob_ana) {
        draw_ana <- sample(1:natt, size = round(runif(n = 1, min = 1, max = natt - 1)), replace = FALSE)
        for (att in 1:natt) {
          if (att %in% draw_ana) mat_ana[resp, ((att * nlevel - att) - 1):(att * nlevel - att)] <- 0
        }
      }
      
      # If screening is flagged, with prob_screen, simulate screening where a respondent screens based 
      # on at least one attribute level but not all of them.
      if (ind_screen == TRUE & runif(1) < prob_screen) {
        draw_screen <- sample(1:nbeta, size = round(runif(n = 1, min = 1, max = nbeta - 1)), replace = FALSE)
        mat_screen[resp, draw_screen] <- -100
      }
      
      # If respondent quality is flagged, with prob_qual, simulate respondent quality where a random
      # value is added to a respondent's latent utility.
      if (ind_qual == TRUE & runif(1) < prob_qual) {
        draw_qual <- sample(1:5, size = 1)
        mat_qual[resp, 1] <- draw_qual
      }
    }
    
    # If heterogeneity isn't flagged, have the first iteration of pathologies apply for all members.
    if (ind_hetero == FALSE) {
      vec_ana <- mat_ana[1,]
      vec_screen <- mat_screen[1,]
      vec_qual <- mat_qual[1,]
      for (resp in 1:nresp) {
        mat_ana[resp,] <- vec_ana
        mat_screen[resp,] <- vec_screen
        mat_qual[resp,] <- vec_qual
      }
    }
    
    # Generate respondent-level betas as a deviation from the population average, gamma, conditioned
    # on the simulated pathologies. Use the respondent-level betas to produce choice data.
    Gamma <- matrix(runif(nbeta, -1, 2), ncol = nbeta)
    Beta <- matrix(double(nbeta * nresp), ncol = nbeta)
    X <- array(double(nresp * ntask * nalt * nbeta), dim = c(nresp, ntask, nalt, nbeta))
    Y <- matrix(double(nresp * ntask), ncol = ntask)
    for (resp in 1:nresp) {
      # Randomly draw a specific version of design_dummy for this respondent along with
      # the associated pathology matrices to modify the utility function.
      resp_ver <- sample(1:nversion, 1)
      resp_design <- design_dummy[which(design_dummy[,1] == resp_ver),]
      resp_ana <- mat_ana[resp,]
      resp_screen <- mat_screen[resp,]
      resp_qual <- mat_qual[resp,]
      
      # Generate respondent-level betas (conditioned on pathologies) and generate choice data.
      Beta[resp,] <- (Gamma + rnorm(length(Gamma), 0, 1)) * resp_ana + resp_screen
      for (task in 1:ntask) {
        # Draw the design matrix for this specific task, compute the latent utility function
        # (conditioned on pathologies) and apply the logit function to produce probabilities 
        # of choosing each alternative. Use the choice probabilities to make a choice.
        task_design <- resp_design[which(resp_design[,2] == task), 4:(ncol(resp_design))]
        resp_util <- as.vector(task_design %*% Beta[resp,]) + resp_qual
        prob_y <- exp(resp_util) / sum(exp(resp_util))
        Y[resp, task] <- sample(1:nalt, 1, prob = prob_y)
        X[resp, task,,] <- matrix(task_design, nrow = nalt, ncol = ncol(task_design))
      }
    }
    
    ####################################################
    # Incorporate clever_randomization.R
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

#############################
# Need to modify pathology matrices for estimation...
# if (ind_test == 1) {
#   mat_ana <- ifelse(ana.mat == 0, 1, 0)
#   mat_screen <- ifelse(screen.mat != 0, 1, 0)
# }
#############################

