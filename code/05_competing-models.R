# Preamble ----------------------------------------------------------------
# Load packages.
library(tidyverse)
library(Rcpp)
library(mvtnorm)
library(bayesm)
library(RcppArmadillo)

# Source functions.
source(here::here("Code", "Source", "ana_hmnl.R"))
source(here::here("Code", "Source", "conj_hmnl_e.R"))

# Set the simulation seed.
set.seed(42)

# Load data and ensemble fit.
data <- read_rds(here::here("Data", str_c("sim_", file_id, ".rds")))
data$train_Z <- matrix(rep(1, nrow(data$train_Y)), ncol = 1)

# Run Competing Models ----------------------------------------------------
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

# Fit attribute non-attendance model.
if (ind_ana == 1) {
  if (!file.exists(here::here("Output", str_c("ana-fit_", file_id, ".rds"))) {
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
    outcomp = ana_hmnl(Data,Prior,Mcmc)
    
    # Save model output.
    write_rds(outcomp, here::here("Output", str_c("ana-fit_", file_id, ".rds")))
  }
}
  
# Fit screening model.
if (ind_screen == 1) {
  if (!here::here("Output", str_c("screen-fit_", file_id, ".rds"))) {
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
    outcomp = conj_hmnl_e(Data,Prior,Mcmc,Cont)
    
    # Save model output.
    write_rds(outcomp, here::here("Output", str_c("screen-fit_", file_id, ".rds")))
  }
}

