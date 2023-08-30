# Preamble ----------------------------------------------------------------
# Load packages.
library(tidyverse)

# Source functions.
source(here::here("code", "source", "simulate_data.R"))
source(here::here("code", "source", "clever_randomization.R"))

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
      hetero = ind_hetero  # Pathologies differ by individual flag.
    )
    
    data <- clever_randomization(
      Y = data$Y,          # Choice data to cleverly randomize.
      X = data$X,          # Design matrices to cleverly randomize.
      natt = 5,            # Number of attributes across design matrices.
      nlevels = rep(3, 5), # Vector of number of attribute levels for each attribute.
      pct_test = .20,      # Percent of data to be saved for testing.
      nmember = 2000       # Number of possible members in the ensemble.
    )
    
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
    
    # Induce clever randomization.
    data <- clever_randomization(
      Y = Y,                                  # Choice data to cleverly randomize.
      X = X,                                  # Design matrices to cleverly randomize.
      natt = 9,                               # Number of attributes across design matrices.
      nlevels = c(5, 8, 3, 4, 2, 2, 2, 2, 2), # Vector of attribute levels for each attribute.
      pct_test = .1,                          # Percent of data to be saved for testing.
      nmember = 2000                          # Number of possible members in the ensemble.
    )
    
    # Save data.
    write_rds(data, here::here("data", str_c("emp_", data_id, ".rds")))
  }
}

