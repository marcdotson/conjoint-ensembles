# Preamble ----------------------------------------------------------------
# This script is used to initialize and run subsequent scripts in the conjoint
# ensembles workflow, including testing parameter and constraint recovery, 
# simulation experiments, and estimating models using real data.

# Load packages.
library(tidyverse)
library(AlgDesign)

# Set the simulation seed.
set.seed(40)

# Use indicator flags to initialize the conjoint ensemble.
ind_none <- 0       # Control for no pathologies.
ind_ana <- 1        # Control for attribute non-attendance.
ind_screen <- 1     # Control for screening.
ind_qual <- 1       # Control for respondent quality.
ind_hetero <- 1     # Control for heterogeneous pathologies.
ind_sim <- 1        # Run a simulation experiment.
ind_test <- 1       # Test for parameter and constraint recovery.
ind_emp <- 0        # Use real empirical data.
ind_beef <- 0       # Use Ground Beef data.
ind_zero <- 0       # Use Zerorez data.

# Specify the arguments for the ensemble members.
nresp <- 300        # Number of respondents.
ntask <- 12         # Number of choice tasks.
nalt <- 3           # Number of choice alternatives.
natt <- 5           # Number of (discrete) attributes.
nlevel <- 3         # Number of attribute levels for each attribute.
nversion <- 10      # Number of versions of the experimental design.
nmember <- 1000     # Number of ensemble members.
prob_ana <- 0.75    # Probability of ANA.
prob_screen <- 0.75 # Probability of screening.
prob_qual <- 0.75   # Probability of respondent quality.
pct_train <- 0.80   # Percent of data used for training.

# Construct the file_id conditioned on flags.
if (ind_none == 1) file_id <- "none"
if (ind_ana == 1) file_id <- "ana"
if (ind_screen == 1) file_id <- "screen"
if (ind_qual == 1) file_id <- "qual"
if (ind_ana == 1 & ind_screen == 1) file_id <- "ana-screen"
if (ind_ana == 1 & ind_qual == 1) file_id <- "ana-qual"
if (ind_screen == 1 & ind_qual == 1) file_id <- "screen-qual"
if (ind_ana == 1 & ind_screen == 1 & ind_qual == 1) file_id <- "ana-screen-qual"
if (ind_hetero == 1) file_id <- paste(file_id, "-hetero", sep = "")
if (ind_hetero == 0) file_id <- paste(file_id, "-homo", sep = "")
if (ind_test == 1) file_id <- paste(file_id, "-test", sep = "")
if (ind_beef == 1) file_id <- paste("ground-beef", "_", file_id, sep = "")
if (ind_zero == 1) file_id <- paste("zerorez", "_", file_id, sep = "")

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

# Upper bounds using the transformed parameters block in hmnl_ensemble_01.
upper_bounds_01 <- tibble(
  Model = rep(c("HMNL", "Ensemble Upper Bound"), 8),
  Pathologies = rep(c(rep("None", 2), rep("ANA", 2), rep("Screen", 2), rep("ANA & Screen", 2)), 2),
  Heterogeneous = c(rep("No", 8), rep("Yes", 8)),
  LOO = NA,
  "Hit Rate" = c(0.582, 0.578, 0.450, 0.457, 0.875, 0.867, 0.911, 0.910,
                 0.582, 0.578, 0.410, 0.421, 0.546, 0.556, 0.594, 0.569),
  "Hit Prob" = c(0.484, 0.482, 0.384, 0.383, 0.862, 0.854, 0.892, 0.881,
                 0.484, 0.482, 0.369, 0.372, 0.525, 0.552, 0.551, 0.566)
)

# Upper bounds using the generated quantities block in hmnl_ensemble_02.
upper_bounds_02 <- tibble(
  Model = rep(c("HMNL", "Ensemble Upper Bound"), 8),
  Pathologies = rep(c(rep("None", 2), rep("ANA", 2), rep("Screen", 2), rep("ANA & Screen", 2)), 2),
  Heterogeneous = c(rep("No", 8), rep("Yes", 8)),
  LOO = NA,
  "Hit Rate" = c(0.582, 0.578, 0.450, 0.450, 0.875, 0.865, 0.911, 0.903,
                 0.582, 0.578, 0.410, 0.410, 0.546, 0.551, 0.594, 0.569),
  "Hit Prob" = c(0.484, 0.482, 0.384, 0.383, 0.862, 0.851, 0.892, 0.883,
                 0.484, 0.482, 0.369, 0.368, 0.525, 0.548, 0.551, 0.565)
)

upper_bounds_01 |> 
  ggplot(aes(x = Model, y = `Hit Prob`)) +
  geom_col() +
  facet_grid(Pathologies ~ Heterogeneous)

upper_bounds_02 |> 
  ggplot(aes(x = Model, y = `Hit Prob`)) +
  geom_col() +
  facet_grid(Pathologies ~ Heterogeneous)

