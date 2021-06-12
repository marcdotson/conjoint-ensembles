# Set Indicator Flags -----------------------------------------------------
# Indicate simulated or empirical data.
ind_sim <- 1 # Indicates simulated data.
ind_emp <- 0 # Indicates empirical data.

if (ind_sim == 1) {
  # Indicate which pathologies to simulate.
  ind_none <- 0       # Indicates no pathologies.
  ind_ana <- 1        # Indicates attribute non-attendance.
  ind_screen <- 0     # Indicates screening.
  ind_resp <- 1       # Indicates respondent quality (bootstrap).
  
  # Decide on pathology heterogeneity and the size of the ensemble.
  ind_hetero <- 1     # Indicates if pathologies differ by individual.
  nmember <- 1000     # Indicates the number of ensemble members.
  
  # Construct the file_id conditioned on flags.
  if (ind_none == 1) file_id <- "none"
  if (ind_ana == 1) file_id <- "ana"
  if (ind_screen == 1) file_id <- "screen"
  if (ind_resp == 1) file_id <- "resp"
  if (ind_ana == 1 & ind_screen == 1) file_id <- "ana-screen"
  if (ind_ana == 1 & ind_resp == 1) file_id <- "ana-resp"
  if (ind_screen == 1 & ind_screen == 1) file_id <- "screen-resp"
  if (ind_ana == 1 & ind_screen == 1 & ind_resp == 1) file_id <- "ana-screen-resp"
  if (ind_hetero == 1) file_id <- paste(file_id, "-hetero", sep = "")
  if (ind_hetero == 0) file_id <- paste(file_id, "-homo", sep = "")
}
if (ind_emp == 1) {
  # Indicate which empirical data to use.
  ind_beef <- 1       # Indicates Ground Beef.
  ind_zero <- 0       # Indicates Zerorez.
  
  # Indicate which pathologies to randomize for.
  # SHOULD BE ALL BY DEFAULT!!!
  ind_ana <- 1        # Indicates attribute non-attendance.
  ind_screen <- 0     # Indicates screening.
  ind_resp <- 1       # Indicates respondent quality (bootstrap).
  
  # Decide on pathology heterogeneity and the size of the ensemble.
  ind_hetero <- 1     # Indicates if pathologies differ by individual.
  nmember <- 1000     # Indicates the number of ensemble members.
  
  # Construct the file_id conditioned on flags.
  if (ind_beef == 1) file_id <- "ground-beef"
  if (ind_zero == 1) file_id <- "zerorez"
  # if (ind_ana == 1) file_id <- paste(file_id, "_ana", sep = "")
  # if (ind_screen == 1) file_id <- paste(file_id, "_screen", sep = "")
  # if (ind_resp == 1) file_id <- paste(file_id, "_resp", sep = "")
  if (ind_ana == 1 & ind_screen == 1) file_id <- paste(file_id, "_ana-screen", sep = "")
  if (ind_ana == 1 & ind_resp == 1) file_id <- paste(file_id, "_ana-resp", sep = "")
  if (ind_screen == 1 & ind_resp == 1) file_id <- paste(file_id, "_screen-resp", sep = "")
  if (ind_ana == 1 & ind_screen == 1 & ind_resp == 1) file_id <- paste(file_id, "_ana-screen-resp", sep = "")
  if (ind_hetero == 1) file_id <- paste(file_id, "-hetero", sep = "")
  if (ind_hetero == 0) file_id <- paste(file_id, "-homo", sep = "")
}

# Run the Ensemble and Competing Models -----------------------------------
# Simulate data or clean empirical data and induce randomization.
source(here::here("Code", "02_data-prep.R"))

# Run the conjoint ensemble using the clever randomization.
source(here::here("Code", "03_conjoint-ensemble.R"))

# Produce weights using the ensemble output.
source(here::here("Code", "04_meta-learner.R"))

# # Run the models specific to the indicated pathology.
# source(here::here("Code", "05_competing-models.R"))

# Compute and compare fit across models.
source(here::here("Code", "06_model-comparison.R"))
model_comparison

