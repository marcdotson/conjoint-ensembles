# Preamble ----------------------------------------------------------------
# Load packages.
library(tidyverse)

# Source functions.
source(here::here("Code", "Source", "simulate_data.R"))
source(here::here("Code", "Source", "clever_randomization.R"))

# # Simulate Data -----------------------------------------------------------
# ind_none <- 0       # Indicates no pathologies.
# ind_ana <- 1        # Indicates attribute non-attendance.
# ind_screen <- 1     # Indicates screening.
# ind_ana_screen <- 0 # Indicates attribute non-attendance and screening.
# 
# hetero <- 1         # Indicates if pathologies differ by individual 
# 
# if (ind_none == 1) file_name <- "none"
# if (ind_ana == 1) file_name <- "ana"
# if (ind_screen == 1) file_name <- "screen"
# if (ind_ana_screen == 1) file_name <- "ana-screen"
# if (hetero == 1) file_name <- paste(file_name,"-hetero", sep="")
# if (hetero == 0) file_name <- paste(file_name,"-homo", sep="")

# Set the simulation seed.
set.seed(42)

ana_logical <- FALSE
screen_logical <- FALSE
hetero_logical <- FALSE
if(ind_ana==1) ana_logical <- TRUE
if(ind_screen==1) screen_logical <- TRUE
if(ind_ana_screen==1) ana_logical <- TRUE
if(ind_ana_screen==1) screen_logical <- TRUE
if(hetero==1) hetero_logical <- TRUE


# Simulate data, impose clever randomization, and save.
data <- simulate_data(
  nhh = 300,      # Number of respondents (households).
  ntask = 12,     # Number of choice tasks.
  nalt = 3,       # Number of choice alternatives.
  natt = 5,       # Number of attributes.
  nlevel = 3,     # Number of estimable attribute levels for each attribute.
  nversion = 10,  # Number of versions of the experimental design.
  ana = ana_logical,    # Attribute non-attendance flag.
  screen = screen_logical, # Screening pathology flag.
  hetero = hetero_logical  # heterogeneous pathologies
)

data <- clever_randomization(
  Y = data$Y,     # Choice data to cleverly randomize.
  X = data$X,     # Design matrices to cleverly randomize.
  pct_test = .20, # Percent of data to be saved for testing.
  nmember = 250 # Number of members in the ensemble.
)

# Save simulated data.
write_rds(data, here::here("Data", str_c("sim_", file_name, ".rds")))

