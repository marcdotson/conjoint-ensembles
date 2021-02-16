## This file controls the run flow for individual scripts 

# Simulate Data -----------------------------------------------------------
ind_none <- 0       # Indicates no pathologies.
ind_ana <- 0        # Indicates attribute non-attendance.
ind_screen <- 1     # Indicates screening.
ind_ana_screen <- 0 # Indicates attribute non-attendance and screening.

hetero <- 1         # Indicates if pathologies differ by individual 

if (ind_none == 1) file_name <- "none"
if (ind_ana == 1) file_name <- "ana"
if (ind_screen == 1) file_name <- "screen"
if (ind_ana_screen == 1) file_name <- "ana-screen"
if (hetero == 1) file_name <- paste(file_name,"-hetero", sep="")
if (hetero == 0) file_name <- paste(file_name,"-homo", sep="")

# Run Ensemble-----------------------------------------------------------

source(here::here("Code","01_simulate-data-hetero.R"))

source(here::here("Code","03_conjoint-ensemble.R"))

source(here::here("Code","04_meta-learner.R"))

source(here::here("Code","06_model-comparison.R"))

model_comparison
