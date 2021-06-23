# Set Indicator Flags -----------------------------------------------------
# Indicate simulated or empirical data.
ind_sim <- 1        # Indicates simulated data.
ind_emp <- 0        # Indicates empirical data.

if (ind_sim == 1) {
  # Indicates a test where we pass the actual constraint 
  # matrices into the ensemble estimation.
  ind_test <- 1
}

if (ind_emp == 1) {
  # Indicate which empirical data to use.
  ind_beef <- 1     # Indicates Ground Beef.
  ind_zero <- 0     # Indicates Zerorez.
}
  
# Indicate which pathologies to randomize for.
ind_none <- 0       # Indicates no pathologies.
ind_ana <- 1        # Indicates attribute non-attendance.
ind_screen <- 0     # Indicates screening.
ind_resp <- 0       # Indicates respondent quality (bootstrap).

# Decide on pathology heterogeneity and the size of the ensemble.
ind_hetero <- 0     # Indicates if pathologies differ by individual.
nmember <- 2     # Indicates the number of ensemble members.

# Construct the file_id conditioned on flags.
if (ind_none == 1) file_id <- "none"
if (ind_ana == 1) file_id <- "ana"
if (ind_screen == 1) file_id <- "screen"
if (ind_resp == 1) file_id <- "resp"
if (ind_ana == 1 & ind_screen == 1) file_id <- "ana-screen"
if (ind_ana == 1 & ind_resp == 1) file_id <- "ana-resp"
if (ind_screen == 1 & ind_resp == 1) file_id <- "screen-resp"
if (ind_ana == 1 & ind_screen == 1 & ind_resp == 1) file_id <- "ana-screen-resp"
if (ind_hetero == 1) file_id <- paste(file_id, "-hetero", sep = "")
if (ind_hetero == 0) file_id <- paste(file_id, "-homo", sep = "")

# Finalize the file_id conditioned on flags.
if (ind_sim == 1) {
  if (ind_test == 1) file_id <- paste(file_id, "-test", sep = "")
}

if (ind_emp == 1) {
  if (ind_beef == 1) data_id <- "ground-beef"
  if (ind_zero == 1) data_id <- "zerorez"
  file_id <- paste(data_id, "_", file_id, sep = "")
}

# Run the Ensemble and Competing Models -----------------------------------
# Simulate data or clean empirical data and induce randomization.
source(here::here("Code", "02_data-prep.R"))

# Run the conjoint ensemble using the clever randomization.
source(here::here("Code", "03_conjoint-ensemble.R"))

# # Produce weights using the ensemble output.
# source(here::here("Code", "04_meta-learner.R"))

# # Run the models specific to the indicated pathology.
# source(here::here("Code", "05_competing-models.R"))

# Compute and compare fit across models.
source(here::here("Code", "06_model-comparison.R"))

ensemble_fit$ensemble_weights
model_comparison

# Parameter recovery...

