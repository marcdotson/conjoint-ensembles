clever_randomization <- function(Y = NULL, X = NULL, pct_test = .20, nmember = 100){
  # This function splits and returns training and testing data and 
  # two matrices that indicates if an attribute level is used for 
  # screening or is not attended to.
  
  # Needs to be modified to all for control over the type of pathology investigated 
  # and the number of attribute levels used for both screening and ana.
  
  # Compute sample characteristics.
  nobs = dim(X)[1] # Number of observational units.
  natt = dim(X)[4] # Number of attribute levels (after dummy coding).
  
  # Split the data.
  ntest = round(pct_test * nobs, 0)
  test_ind = sort(sample(1:nobs, ntest, replace = FALSE))
  '%!in%' <- function(x, y)!('%in%'(x, y))
  train_ind = which(c(1:nobs) %!in% test_ind)
  
  train_Y = Y[train_ind,]
  train_X = X[train_ind,,,]
  
  test_Y = Y[test_ind,]
  test_X = X[test_ind,,,]
  
  # Create randomization patterns for Screening and ANA
  # Allow for 1 screened attribute level and 3 ANA levels
  # Repeat for the number of ensemble members
  mat_ana = mat_screen = matrix(double(nmember * natt), ncol = natt)
  for(i in 1:nmember){
    ana.ind = sample(1:natt,3,replace = FALSE)
    screen.ind = sample(1:natt,1,replace = FALSE)
    mat_ana[i,ana.ind] = 1
    mat_screen[i,screen.ind] = 1
  }
  
  return(list(train_Y = train_Y, train_X = train_X, test_Y = test_Y,
              test_X = test_X, mat_ana = mat_ana, mat_screen = mat_screen))
}
