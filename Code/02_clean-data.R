## Clean data for use in the empirical application. 

## This code should live somewhere else.  Including it here so I can test the clean_mnl funcion locally
## Using real datasets for testing.  Need to modifiy the simulation routines to generate data with the same structure.
## Using zerorez data from my conjoint datasets file
## This needs to be written as a function

## Read Design Data
desmat = read.csv(here::here("Data","coded.design.csv"),header=TRUE)
## Read Choice Data
choice.dat = read.csv(here::here("Data","zerorez.final.csv"),header=TRUE)
choice.dat = choice.dat[,-1]

## Sample Characteristics
nobs = nrow(choice.dat)
ntask = 10
nalt = 4
natt = ncol(desmat) - 3

## Format Data for Stan
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


## This function splits the data into training and test subsets
## It returns 2 matricies of dimnesion n.member * natt that indicate if 
## an attribute level is used for screening or is not attended to
## It also returns the training and test datasets
## Needs to be modified to all for control over the type of pathology investigated 
## and the number of attribute levels used for both screening and ana

clean_mnl <- function(Y=NULL, X=NULL, pct.test = .20, n.member = 100){

  ## compute sample characteristics
  nobs = dim(X)[1] #number of observational units
  natt = dim(X)[4] #number of attribute levels (after dummy coding)
  
  ## Split the data
  n.test = round(pct.test * nobs,0)
  ind.test = sort(sample(1:nobs,n.test,replace =FALSE))
  '%!in%' <- function(x,y)!('%in%'(x,y))
  ind.train = which(c(1:nobs) %!in% ind.test)
  
  Y.train = Y[ind.train,]
  X.train = X[ind.train,,,]
  
  Y.test = Y[ind.test,]
  X.test = X[ind.test,,,]
  
  ## Create randomization patterns for Screening and ANA
  ## Allow for 1 screended attribute level and 3 ANA levels
  ## Repeat for the number of ensemble members
  
  ana.mat = screen.mat = matrix(double(n.member*natt),ncol = natt)
  
  for(i in 1:n.member){
    ana.ind = sample(1:natt,3,replace = FALSE)
    screen.ind = sample(1:natt,1,replace = FALSE)
    ana.mat[i,ana.ind] = 1
    screen.mat[i,screen.ind] = 1
  }
  
  return(list(Y.train = Y.train, X.train = X.train, Y.test = Y.test,
              X.test = X.test, ana.mat = ana.mat, screen.mat = screen.mat))
}


## Run clean_mnl code

XX = clean_mnl(Y = Y, X = X, pct.test = .1, n.member = 200)

save(XX, file = here::here("Data", "design.RData"))

