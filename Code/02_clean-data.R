# This code should live somewhere else.  Including it here so I can test the
# clever_randomization function locally Using real datasets for testing.  Need to modify
# the simulation routines to generate data with the same structure. Using
# zerorez data from my conjoint datasets file

# This needs to be written as a function

# Read Design Data
#desmat <- read.csv(here::here("Data", "R05_Zerorez", "coded.design.csv"),header=TRUE)
desmat <- read.csv(here::here("Data", "R06_Ground-Beef", "coded.design.gb.csv"),header=TRUE)
desmat <- desmat[,-1]
# Read Choice Data
#choice.dat = read.csv(here::here("Data", "R05_Zerorez", "zerorez.final.csv"),header=TRUE)
choice.dat = read.csv(here::here("Data", "R06_Ground-Beef", "jcb.results1.csv"),header=TRUE)
choice.dat = choice.dat[,-1]

# Sample Characteristics
nobs = nrow(choice.dat)
ntask = 10
nalt = 4
natt = ncol(desmat) - 3

# Format Data for Stan
regdata = NULL
X = array(double(nobs*ntask*nalt*natt),dim=c(nobs,ntask,nalt,natt))
Y = as.matrix(choice.dat[,c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10")])

for(i in 1:nobs){
  #X = array(double(ntask*nalt*natt),dim=c(nalt,natt,ntask))
  #Y = as.vector(choice.dat[i,c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10")])
  #Y = as.vector(as.matrix(Y))
  nver = choice.dat[i,"X0.Version"]
  tmp = desmat[which(desmat[,1]==nver),]
  for(j in 1:ntask){
    X[i,j,,] = as.matrix(tmp[which(tmp[,2]==j),4:(ncol(desmat))])
  }
}

# Run clever_randomization code
source(here::here("Code", "Source", "clever_randomization.R"))
data = clever_randomization(Y = Y, X = X, pct_test = .1, nmember = 10)

save(data, file = here::here("Data", "R06_Ground-Beef", "design.RData"))

