# Preamble ----------------------------------------------------------------
# This script is used to initialize and run subsequent scripts in the conjoint
# ensembles workflow, including testing parameter and constraint recovery, 
# simulation experiments, and estimating models using real data.

# Load packages.
library(AlgDesign)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(tidybayes)
library(loo)

# Set the simulation seed, loo option.
set.seed(40)
options(mc.cores = 4)

# Use indicator flags to initialize the conjoint ensemble.
ind_test <- 1       # Test for parameter and constraint recovery.
ind_sim <- 1        # Run a simulation experiment.
# ind_emp <- 0        # Use real empirical data.
# ind_beef <- 0       # Use Ground Beef data.
# ind_zero <- 0       # Use Zerorez data.
ind_none <- 0       # Control for no pathologies.
ind_ana <- 1        # Control for attribute non-attendance.
ind_screen <- 1     # Control for screening.
ind_qual <- 1       # Control for respondent quality.
ind_hetero <- 1     # Control for heterogeneous pathologies.

# Specify the arguments for the ensemble members.
nresp <- 300        # Number of respondents.
ntask <- 12         # Number of choice tasks.
nalt <- 3           # Number of choice alternatives.
natt <- 5           # Number of (discrete) attributes.
nlevel <- 3         # Number of attribute levels for each attribute.
ncov <- 0           # Number of population-level covariates.
nversion <- 10      # Number of versions of the experimental design.
nmember <- 1        # Number of ensemble members.
prob_ana <- 0.50    # Probability of ANA.
prob_screen <- 0.75 # Probability of screening.
prob_qual <- 0.20   # Probability of respondent quality.
pct_train <- 0.80   # Percent of data used for training.

# Construct the patho_id conditioned on flags.
if (ind_sim == 1) data_id <- "sim"
# if (ind_emp == 1) data_id <- "emp"
# if (ind_beef == 1) data_id <- str_c(data_id, "_", "ground-beef")
# if (ind_zero == 1) data_id <- str_c(data_id, "_", "zerorez")
if (ind_none == 1) patho_id <- "none"
if (ind_ana == 1) patho_id <- "ana"
if (ind_screen == 1) patho_id <- "screen"
if (ind_qual == 1) patho_id <- "qual"
if (ind_ana == 1 & ind_screen == 1) patho_id <- "ana-screen"
if (ind_ana == 1 & ind_qual == 1) patho_id <- "ana-qual"
if (ind_screen == 1 & ind_qual == 1) patho_id <- "screen-qual"
if (ind_ana == 1 & ind_screen == 1 & ind_qual == 1) patho_id <- "ana-screen-qual"
if (ind_hetero == 0) patho_id <- str_c(patho_id, "-homo")
if (ind_hetero == 1) patho_id <- str_c(patho_id, "-hetero")
if (ind_test == 1) patho_id <- str_c(patho_id, "_test")

data_id
patho_id

# Run the Ensemble and Competing Models -----------------------------------
# Prepare choice data and induce clever randomization.
source(here::here("code", "02_data-prep.R"))

# Run the conjoint ensemble using clever randomization.
source(here::here("code", "03_conjoint-ensemble.R"))

# # Produce weights using the ensemble output.
# source(here::here("code", "04_meta-learner.R"))

# # Run the models specific to the indicated pathology.
# source(here::here("code", "05_competing-models.R"))

# Compute and compare fit across models.
source(here::here("code", "06_model-comparison.R"))
# ensemble_fit$ensemble_weights
model_comparison
