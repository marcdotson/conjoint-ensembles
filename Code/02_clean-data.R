# Preamble ----------------------------------------------------------------
# Load packages.
library(tidyverse)

# Source functions.
source(here::here("Code", "Source", "clever_randomization.R"))

# Load Data ---------------------------------------------------------------
ind_beef <- 1       # Indicates Ground Beef.
ind_zero <- 0       # Indicates Zerorez.

if (ind_beef == 1) file_name <- "ground-beef"
if (ind_zero == 1) file_name <- "zerorez"

# Load design.
design <- read.csv(here::here("Data", str_c("emp_", file_name, "_design.csv")), header = TRUE)
design <- design[,-1]

# Load choice data.
choice_data <- read.csv(here::here("Data", str_c("emp_", file_name, "_final.csv")), header=TRUE)
choice_data <- choice_data[,-1]

# Clean Data and Induce Randomization -------------------------------------
# Sample Characteristics
nobs = nrow(choice_data)
ntask = 10
nalt = 4
natt = ncol(design) - 3

# Format Data for Stan
regdata = NULL
X = array(double(nobs*ntask*nalt*natt),dim=c(nobs,ntask,nalt,natt))
Y = as.matrix(choice_data[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10")])

for(i in 1:nobs){
  nver = choice_data[i,"X0.Version"]
  tmp = design[which(design[,1]==nver),]
  for(j in 1:ntask){
    X[i,j,,] = as.matrix(tmp[which(tmp[,2]==j),4:(ncol(design))])
  }
}

# Induce randomization.
data = clever_randomization(Y = Y, X = X, pct_test = .1, nmember = 10)

# Save data.
write_rds(data, here::here("Data", str_c("emp_", file_name, ".rds")))

