# Preamble ----------------------------------------------------------------
# Load packages.
library(tidyverse)
library(Rcpp)
library(mvtnorm)
library(bayesm)
library(RcppArmadillo)

# Load Data ---------------------------------------------------------------
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

data <- read_rds(here::here("Data", str_c("sim_", file_name, ".rds")))
Y <- data$train_Y
X <- data$train_X
Z <- matrix(rep(1, nrow(Y)), ncol = 1)


# ...the rest needs to be updated.
# QUESTIONS:
# - What functions are we using for the competing models?
# - Update saving output using file_name.


# Array of experimental designs per choice task.
#X_raw <- read_csv(here::here("Data", "05_RealDate", "X.csv"))

nvar <- (ncol(X_raw)-3)
nz <- ncol(Z)
nresp <- nrow(Y)
nalt <- nrow(X_raw)/length(Y)
ntask <- nrow(X_raw)/(nalt*nresp)

# Attribute Levels
nlvl <- rep(2,nvar)                            # Number of levels for each attribute (assuming all binary here)
lvlstart <- cumsum(nlvl)-(nlvl-1) - c(0:(length(nlvl)-1))

# MCMC parameters
R <- 1000                                              # Number of iterations in the Markov chain
keep <- 5                                              # Thinning parameter
print <- 100                                           # Skip interval for printing output
step <- 2.93/(nvar^.5)                                 # RW step size


# Priors
gammabar <- matrix(rep(0,nz*nvar),nrow=nz)              # Vector of means for Gamma
Agamma <- 0.01*diag(nz)                                 # Gamma precision matrix
V <- (nvar+3)*diag(nvar)                                # Location for IW prior for Vbeta
nu <- nvar+3                                            # DoF for IW prior for Vbeta


# Fit screening model
if(m_ind_scr==1){
  # Create input lists
  data=NULL
  for(resp in 1:nresp){
    y <- Y[resp,] %>% 
      as.matrix()
    X <- X_raw[((resp-1)*nalt*ntask + 1):(resp*nalt*ntask),4:ncol(X_raw)]%>%
      as.matrix()
    S <- NULL
    for(att in 1:length(nlvl)){
      X_att <- X[,lvlstart[att]:cumsum(nlvl-1)[att]]%>%
        as.matrix()
      S_att <- cbind(rep(1, nrow(X)) - rowSums(X_att),X_att)
      S <- cbind(S,S_att)
    }
    data[[resp]]=list(y=y, X=X, S=S )
  }
  
  Data <- list(nalts=nalt, data=data, W=Z, Beta=Beta)
  Prior <- list(gammabar=gammabar, Agamma=Agamma, nu=nu, V=V)
  Mcmc <- list(R=R, keep=keep, print=print, bstep=step, sim_ind=0, cont_ind=ind_cont, het_ind=0, out_ind=ind_out)
  Cont=NULL
  
  # Fit screening model to data
  source("Code/R/Source/conj_hmnl_e.R"); outcomp = conj_hmnl_e(Data,Prior,Mcmc,Cont)
  
  if(ind_non==1) write_rds(outcomp, here::here("Output", "Competing Models", "conj_hmnl_on_nonData.RDS"))
  if(ind_scr==1) write_rds(outcomp, here::here("Output", "Competing Models", "conj_hmnl_on_scrData.RDS"))
  if(ind_ana==1) write_rds(outcomp, here::here("Output", "Competing Models", "conj_hmnl_on_anaData.RDS"))
  if(ind_sna==1) write_rds(outcomp, here::here("Output", "Competing Models", "conj_hmnl_on_snaData.RDS"))
}

# Fit attribute non-attendance model
if(m_ind_ana==1){
  # Create input lists
  data=NULL
  for(resp in 1:nresp){
    data[[resp]]=list(y=Y[resp,] %>% 
                        as.matrix(),
                      X_raw=X_raw[((resp-1)*nalt*ntask + 1):(resp*nalt*ntask),]%>%
                        as.matrix(),
                      S=S[resp,] %>%
                        as.matrix())
  }
  alpha=c(2,2)
  Data <- list(nalts=nalt, data=data, Z=Z, lvlvec=nlvl)
  Prior <- list(Deltabar=gammabar, ADelta=Agamma, nu=nu, V=V, alpha=alpha)
  Mcmc <- list(R=R, keep=keep, print=print, bstep=bstep, sim_ind=0, cont_ind=ind_cont, het_ind=0)
  
  # Fit screening model to data
  source("Code/R/Source/ana_hmnl.R"); outcomp = ana_hmnl(Data,Prior,Mcmc)
  
  if(ind_non==1) write_rds(outcomp, here::here("Output", "Competing Models", "ana_hmnl_on_nonData.RDS"))
  if(ind_scr==1) write_rds(outcomp, here::here("Output", "Competing Models", "ana_hmnl_on_scrData.RDS"))
  if(ind_ana==1) write_rds(outcomp, here::here("Output", "Competing Models", "ana_hmnl_on_anaData.RDS"))
  if(ind_sna==1) write_rds(outcomp, here::here("Output", "Competing Models", "ana_hmnl_on_snaData.RDS"))
}