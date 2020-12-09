# Load packages.
library(tidyverse)
# library(rstan)
# library(bayesplot)
# library(tidybayes)
library(loo)

# Source fit functions.
source(here::here("Code", "Source", "predictive_fit_ensemble.R"))
source(here::here("Code", "Source", "hit_rate.R"))
source(here::here("Code", "Source", "hit_prob.R"))

ind_none <- 0       # Indicates no pathologies.
ind_ana <- 1        # Indicates attribute non-attendance.
ind_screen <- 0     # Indicates screening.
ind_ana_screen <- 0 # Indicates attribute non-attendance and screening.
ind_real <- 0       # Indicates ____ data.

if (ind_none == 1) file_name <- "none"
if (ind_ana == 1) file_name <- "ana"
if (ind_screen == 1) file_name <- "screen"
if (ind_ana_screen == 1) file_name <- "ana-screen"
if (ind_real == 1) file_name <- "design"

# Load hmnl-fit output.
hmnl_fit <- read_rds(here::here("Output", str_c("hmnl-fit_", file_name, ".rds")))

# Compute fit.
loo(hmnl_fit)

# Save as a data frame with ensemble-fit and competing model fit.

# Start with loo.
# Hit rate and hit prob in Code/Source.
