# Set Indicator Flags -----------------------------------------------------
# Indicate simulated or empirical data.
ind_sim <- 1        # Indicates simulated data.
ind_emp <- 0        # Indicates empirical data.

if (ind_sim == 1) {
  # Indicates a test where we pass on the actual constraint 
  # matrices and estimate a single, best-performing model.
  ind_test <- 1
}

if (ind_emp == 1) {
  # Indicate which empirical data to use.
  ind_beef <- 1     # Indicates Ground Beef.
  ind_zero <- 0     # Indicates Zerorez.
}

# Indicate which pathologies to randomize for.
ind_none <- 1       # Indicates no pathologies.
ind_ana <- 0        # Indicates attribute non-attendance.
ind_screen <- 0     # Indicates screening.
ind_resp <- 0       # Indicates respondent quality (still a bootstrap).

# Decide on pathology heterogeneity and the size of the ensemble.
ind_hetero <- 1     # Indicates if pathologies differ by individual.
nmember <- 1000     # Indicates the number of ensemble members.
if (ind_test == 1) {
  nmember <- 1      # Constrain tests to a single, best-performing model.
}

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
if (ind_test == 1) file_id <- paste(file_id, "-test", sep = "")

if (ind_emp == 1) {
  if (ind_beef == 1) data_id <- "ground-beef"
  if (ind_zero == 1) data_id <- "zerorez"
  file_id <- paste(data_id, "_", file_id, sep = "")
}

file_id

# Run the Ensemble and Competing Models -----------------------------------
# Simulate data or clean empirical data and induce randomization.
source(here::here("code", "02_data-prep.R"))

# Run the conjoint ensemble using the clever randomization.
source(here::here("code", "03_conjoint-ensemble.R"))

# # Produce weights using the ensemble output.
# source(here::here("code", "04_meta-learner.R"))

# # Run the models specific to the indicated pathology.
# source(here::here("code", "05_competing-models.R"))

# Compute and compare fit across models.
source(here::here("code", "06_model-comparison.R"))
# ensemble_fit$ensemble_weights
# model_comparison

# Print results.
model_comparison

# Upper bounds.
upper_bounds <- tibble(
  Model = rep(c("HMNL", "Ensemble"), 4),
  Pathologies = c(rep("None", 2), rep("ANA", 2), rep("Screen", 2), rep("ANA & Screen", 2)),
  Heterogeneous = rep("No", 8),
  LOO = NA,
  "Hit Rate" = c(0.578, 0.574, 0.614, 0.604, 0.739, 0.742, 0.865, 0.717),
  "Hit Prob" = c(0.483, 0.482, 0.532, 0.527, 0.685, 0.609, 0.832, 0.682)
)

model_comparison
upper_bounds

