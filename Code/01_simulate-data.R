# Load readr.
library(readr)

# Source the simulate_data() and clever_randomization() functions.
source(here::here("Code", "Source", "simulate_data.R"))
source(here::here("Code", "Source", "clever_randomization.R"))

# Simulate data, impose clever randomization, and save.
data <- simulate_data(
  nhh = 300,     # Number of respondents (households).
  ntask = 12,    # Number of choice tasks.
  nalt = 3,      # Number of choice alternatives.
  natt = 5,      # Number of attributes.
  nlevel = 3,    # Number of estimable attribute levels for each attribute.
  nversion = 10, # Number of versions of the experimental design.
  ana = TRUE,   # Attribute non-attendance flag.
  screen = TRUE # Screening pathology flag.
)

data <- clever_randomization(
  Y = data$Y,     # Choice data to cleverly randomize.
  X = data$X,     # Design matrices to cleverly randomize.
  pct_test = .20, # Percent of data to be saved for testing.
  nmember = 100   # Number of members in the ensemble.
)

write_rds(data, here::here("Data", "sim_ana-screen.rds"))

