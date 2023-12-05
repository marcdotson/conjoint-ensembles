# Prepare choice data and induce clever randomization.

##########################
# Integrate cleaning and reformatting specific empirical data and overwriting simulated
# data to avoid duplication and result in writing the same data structure.
##########################

# if (ind_emp == 1) {
#   # Load design.
#   design <- read_csv(here::here("data", str_c("emp_", data_id, "_design.csv")))
#   design <- design[,-1]
# 
#   # Load choice data.
#   choice_data <- read_csv(here::here("data", str_c("emp_", data_id, "_final.csv")))
#   choice_data <- choice_data[,-1]
# 
#   # Overwrite ensemble arguments with data.
#   nresp <- nrow(choice_data)
#   ntask <- max(design$Task)
#   nalt <- max(design$Concept)
#   nbeta <- ncol(design) - 3
# 
#   # Reformat data for Stan.
#   X <- array(double(nresp * ntask * nalt * nbeta), dim = c(nresp, ntask, nalt, nbeta))
#   Y <- as.matrix(choice_data[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10")])
#   for (i in 1:nobs) {
#     nver <- choice_data[i,"X0.Version"]
#     tmp <- design[which(design[,1]==nver),]
#     for (j in 1:ntask) {
#       X[i,j,,] <- as.matrix(tmp[which(tmp[,2] == j), 4:(ncol(design))])
#     }
#   }
# }

if (!file.exists(here::here("data", str_c(data_id, "_", patho_id, ".rds")))) {
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
  mat_ana <- matrix(double(nresp * nbeta), ncol = nbeta) + 1
  mat_screen <- matrix(double(nresp * nbeta), ncol = nbeta)
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
    
    # If respondent quality is flagged, with prob_qual, simulate respondent quality where a (random?)
    # value is added to a respondent's latent utility.
    if (ind_qual == TRUE & runif(1) < prob_qual) mat_qual[resp, 1] <- 5
  }
  
  # If heterogeneity isn't flagged, have the first iteration of pathologies apply for all respondents.
  if (ind_hetero == 0) {
    tmp_ana <- mat_ana[1,]
    tmp_screen <- mat_screen[1,]
    tmp_qual <- mat_qual[1,]
    for (resp in 1:nresp) {
      mat_ana[resp,] <- tmp_ana
      mat_screen[resp,] <- tmp_screen
      mat_qual[resp,] <- tmp_qual
    }
  }
  
  # Generate respondent-level betas as a deviation from the population average, Gamma, conditioned
  # on the simulated pathologies. Use the respondent-level betas to produce choice data.
  Gamma <- matrix(runif(nbeta, -1, 2), ncol = nbeta)
  Beta <- matrix(double(nresp * nbeta), ncol = nbeta)
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
  
  # Randomly assign respondents into training and testing data.
  '%!in%' <- function(x, y)!('%in%'(x, y))
  nresp_train <- round(pct_train * nresp, 0)
  nresp_test <- round((1 - pct_train) * nresp, 0)
  index_train <- sort(sample(1:nresp, nresp_train, replace = FALSE))
  index_test <- which(c(1:nresp) %!in% index_train)
  
  # Split the choice data and design matrices into training and testing.
  train_Y <- Y[index_train,]
  train_X <- X[index_train,,,]
  test_Y <- Y[index_test,]
  test_X <- X[index_test,,,]
  
  # Generate an array of clever randomization patterns for each possible pathology
  # for each respondent in the training data for each possible member of the ensemble.
  array_ana <- array(double(nbeta * nresp_train * nmember), dim = c(nresp_train, nbeta, nmember))
  array_screen <- array(double(nbeta * nresp_train * nmember), dim = c(nresp_train, nbeta, nmember))
  array_qual <- array(double(nresp_train * nmember), dim = c(nresp_train, 1, nmember))
  for (member in 1:nmember) {
    for (resp_train in 1:nresp_train) {
      # With prob_ana, randomize ANA where a respondent pays attention to at least one attribute
      # and has non-attendance for at least one attribute.
      if (runif(1) < prob_ana) {
        draw_ana <- sample(1:natt, size = round(runif(n = 1, min = 1, max = natt - 1)), replace = FALSE)
        for (att in 1:natt) {
          if (att %in% draw_ana) {
            array_ana[resp_train, ((att * nlevel - att) - 1):(att * nlevel - att), member] <- 1
          }
        }
      }
      
      # With prob_screen, randomize screening where a respondent screens based on at least one attribute
      # level but not all of them.
      if (runif(1) < prob_screen) {
        draw_screen <- sample(1:nbeta, size = round(runif(n = 1, min = 1, max = nbeta - 1)), replace = FALSE)
        array_screen[resp_train, draw_screen, member] <- 1
      }
      
      # With prob_qual, randomize respondent quality where a random value is added to a 
      # respondent's latent utility.
      if (runif(1) < prob_qual) array_qual[resp_train, 1, member] <- 1
    }
  }
  
  # If heterogeneity isn't flagged, have the first iteration of pathologies apply for all respondents
  # in the training data for each member of the pathology.
  if (ind_hetero == 0) {
    for (member in 1:nmember) {
      tmp_ana <- array_ana[1,,member]
      tmp_screen <- array_screen[1,,member]
      tmp_qual <- array_qual[1,,member]
      for (resp_train in 1:nresp_train) {
        array_ana[resp_train,,member] <- tmp_ana
        array_screen[resp_train,,member] <- tmp_screen
        array_qual[resp_train,,member] <- tmp_qual
      }
    }
  }
  
  # If testing is flagged, modify the known simulated pathology matrices into an array
  # with flags for estimation.
  if (ind_test == 1) {
    mat_ana <- ifelse(mat_ana == 0, 1, 0)
    mat_screen <- ifelse(mat_screen != 0, 1, 0)
    mat_qual <- ifelse(mat_qual != 0, 1, 0)
    for (member in 1:nmember) {
      array_ana[,,member] <- mat_ana[index_train,]
      array_screen[,,member] <- mat_screen[index_train,]
      array_qual[,,member] <- mat_qual[index_train,]
    }
  }
  
  # Save the simulated data.
  data <- list(
    # Training and testing data.
    train_Y = train_Y, train_X = train_X, test_Y = test_Y, test_X = test_X,
    # Pathology arrays.
    array_ana = array_ana, array_screen = array_screen, array_qual = array_qual,
    # Population mean Gamma and Beta matrix.
    Gamma = Gamma, Beta = Beta,
    # Index for training and testing.
    index_train = index_train, index_test = index_test
  )
  write_rds(data, here::here("data", str_c(data_id, "_", patho_id, ".rds")))
}

