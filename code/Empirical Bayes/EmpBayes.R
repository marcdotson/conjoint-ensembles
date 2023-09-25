# Create an empty list
emp_data <- vector("list", length = 240)

# Load the necessary R files
source(here::here("code", "01_control-file.R"))
source(here::here("code", "02_data-prep.R"))

# Convert stan data into a list format for empirical data
for (i in 1:240) {
  emp_data[[i]]$y = stan_data$Y[i,]
  emp_data[[i]]$X = stan_data$X[i,,,]
}

## Estimate Empirical Bayes HMNL
rempBayes = function(lgtdata = NULL, lambda = .25){
  
  ###### TODO #####
  ## allow for constraints via lower and upper arguments
  ## more precision, longer tuning for convergence?
  
  require(bayesm)
  #library(mvtnorm)
  
  # This first function calculates the overall log likelihood for different attributes.
  
  ## function to evaluate an MNL with an discrete y
  llmnl.agg = function(beta,lgtdata){
    llike = 0
    J = 10
    for(i in 1:length(lgtdata)){
      X = lgtdata[[i]]$X
      y = lgtdata[[i]]$y
      for(j in 1:J){
        # Calculate the utilities
        u = X[,,j]%*%beta
        # Convert utilities to probabilities
        cprob = exp(u)/sum(exp(u))
        # Take the log og the probablity to update the log likelihood
        llike = llike + log(cprob[y[j]])
      }
    }
    return(llike)
  }
  
  ## function to evaluate an MNL with an allocative y
  
  # This function calculates the individual log likelihood.
  # It blends in the aggregate information.
  llmnl.indv = function(beta,X, y, bbar, lambda = .25){
    llike = 0
    J = 10
    for(j in 1:J){
      # Compute aggregate choice probabilities
      u.agg = X[,,j]%*%bbar
      cprob.agg = as.vector(exp(u.agg)/sum(exp(u.agg)))
      # Create a combined choice vector
      ytmp = double(nrow(X))
      ytmp[y[j]] = 1
      y.comb = ytmp + lambda * cprob.agg
      y.comb = y.comb / sum(y.comb)
      # Compute individual log likelihood
      u.indv = X[,,j]%*%beta
      cprob.indv = as.vector(exp(u.indv) / sum(exp(u.indv)))
      # Take the dot product of the individual probabilities and the combined choice vector
      llike = llike + log(cprob.indv) %*% y.comb
    }
    return(llike)
  }
  
  ## optimize LL for y and X stacked -- Aggregate-level parameters
  # Create empty beta vector
  betainit = double(ncol(lgtdata[[1]]$X[,,1]))
  # Optimize
  out = optim(betainit, llmnl.agg, method = "BFGS", hessian = TRUE, control = list(fnscale = -1, 
                                                                                   trace = 1, reltol = 1e-06), lgtdata = lgtdata)
  H = -1*out$hessian
  #H = mnlHess(out$par,ypooled,Xpooled)
  vbeta = solve(H)
  #vbeta = 100*diag(ncol(H))
  bbar.est = out$par
  
  ## optimize individual level likelihoods and take average (e.g., .25 aggregate)
  betainit = bbar.est 
  
  bmat.est = matrix(double(length(betainit)*length(lgtdata)),ncol=length(betainit))
  
  #lambda = 30
  
  for(i in 1:length(lgtdata)){
    X = lgtdata[[i]]$X
    y = lgtdata[[i]]$y
    # Optimize the individual likelihood using bfgs. Estimate parameters that are individual specific
    # But these parameters are informed by individual and aggregate data.
    outi = optim(betainit, llmnl.indv, method = "BFGS", control = list(fnscale = -1, 
                                                                       trace = 0, reltol = 1e-06), X = X, y = y, bbar=bbar.est, lambda = lambda)
    bmat.est[i,] = outi$par
    #print(i)
  }
  
  # Return the overall output
  return(list(bmat.est = bmat.est, bbar.est = bbar.est))
  
}

# Print the results
rempBayes(emp_data)


