# Compute and compare fit across models.
data <- read_rds(here::here("data", str_c(data_id, "_", patho_id, ".rds")))

####################################################
# We use the average of Beta as the population mean without reference to Sigma.
# Do we just need to import the hmnl_draws rather than hmnl_fit? What about log_lik?
####################################################

hmnl_fit <- read_rds(here::here("output", str_c("hmnl-fit_", data_id, "_", patho_id, ".rds")))
ensemble_draws <- read_rds(here::here("output", str_c("ensemble-fit_", data_id, "_", patho_id, "_", nmember, ".rds")))

# Compute population mean for the choice model with unconstrained betas.
hmnl_draws <- hmnl_fit$draws(format = "df", variables = c("Beta", "Gamma", "Omega", "tau"))
Gamma_mean_01 <- hmnl_draws |> 
  spread_draws(Beta[, j]) |> 
  summarize(mean = mean(Beta)) |> 
  select("mean") |> 
  as.matrix()

# Compute population mean for the choice model with betas constrained after estimation.

####################################################
# Do we need to impose the constraints for each iteration and chain instead?
# Revisit specifying actual pre-built data structure for beta_draws_new
# that accounts for many ensemble members, if needed?
####################################################

beta_draws <- hmnl_draws |> 
  subset_draws(variable = "Beta") |> 
  summarize_draws("mean")
  # as_draws_array()

# beta_draws_new <- array(NA, dim = dim(data$array_ana))
beta_draws_new <- NULL
# for (iter in 1:dim(beta_draws)[1]) {
#   for (chain in 1:dim(beta_draws)[2]) {
    # Restructure the specific draw to match the pathology arrays.
    beta_draws_temp <- matrix(
      # beta_draws[iter, chain,], 
      beta_draws$mean,
      nrow = dim(data$train_X)[1], 
      ncol = dim(data$train_X)[4], 
      byrow = TRUE
    )
    for (member in 1:dim(data$array_ana)[3]) {
      for (resp_train in 1:dim(data$train_X)[1]) {
        # Match the dimensions of the arrays to impose the same sort of constraints
        # on the parameter estimates directly that we do as part of the transformed
        # parameters in the actual ensemble setup.
        for (level in 1:dim(data$train_X)[4]) {
          # Impose fixed values using screening indicator array.
          if (data$array_screen[resp_train,level,member] == 1) {
            beta_draws_temp[resp_train,level] = -100
          }
          
          # Impose fixed values using ANA indicator array.
          if (data$array_ana[resp_train,level,member] == 1) {
            beta_draws_temp[resp_train,level] = 0
          }
        }
        
        # Impose fixed values using respondent quality indicator array.
        if (data$array_qual[resp_train,1,member] == 1) {
          beta_draws_temp[resp_train,] = beta_draws_temp[resp_train,] * 0
        }
        
        # Save modified beta_draws_temp.
        beta_draws_new <- rbind(beta_draws_new, beta_draws_temp[resp_train,])
      }
    }
#   }
# }

Gamma_mean_02 <- matrix(apply(beta_draws_new, 2, mean), ncol = 1)

# Compute population mean for the choice model with betas constrained during estimation.
Gamma_mean_03 <- ensemble_draws[[1]][[1]] |> 
  spread_draws(Beta[, j]) |> 
  summarize(mean = mean(Beta)) |> 
  select("mean") |> 
  as.matrix()

# Create an empty model comparison data frame.
model_comparison <- tibble(
  Model = c("Unconstrained Betas", "Constrained After Estimation", "Constrained During Estimation"),
  ID = str_c(data_id, "_", patho_id),
  LOO = NA,
  "Hit Rate" = NA,
  "Hit Prob" = NA
)

####################
# Problem with loo working with a cmdstanr object?
####################

nresp_test <- dim(data$test_X)[1]
hits_01 <- hits_02 <- hits_03 <- matrix(NA, nrow = nresp_test, ncol = ntask)
probs_01 <- probs_02 <- probs_03 <- matrix(NA, nrow = nresp_test, ncol = ntask)
for (resp_test in 1:nresp_test) {
  for (task in 1:ntask) {
    # Pull the relevant choice y and design matrix X.
    tmp_y <- data$test_Y[resp_test, task]
    tmp_X <- data$test_X[resp_test, task,,]
    
    ####################################
    # In-sample predictive fit would use betas.
    # Will need to use $Beta along with $Gamma.
    ####################################
    
    # Compute model fit for the choice model with unconstrained betas.
    hits_01[resp_test, task] <- tmp_y == which.max(tmp_X %*% Gamma_mean_01)
    tmp_XGamma <- matrix(exp(tmp_X %*% Gamma_mean_01), byrow = TRUE, ncol = nalt)
    probs_01[resp_test, task] <- exp(sum(log(tmp_XGamma) - log(sum(tmp_XGamma))))
    
    # Compute model fit for the choice model with betas constrained after estimation.
    hits_02[resp_test, task] <- tmp_y == which.max(tmp_X %*% Gamma_mean_02)
    tmp_XGamma <- matrix(exp(tmp_X %*% Gamma_mean_02), byrow = TRUE, ncol = nalt)
    probs_02[resp_test, task] <- exp(sum(log(tmp_XGamma) - log(sum(tmp_XGamma))))
    
    # Compute model fit for the choice model with betas constrained during estimation.
    hits_03[resp_test, task] <- tmp_y == which.max(tmp_X %*% Gamma_mean_03)
    tmp_XGamma <- matrix(exp(tmp_X %*% Gamma_mean_03), byrow = TRUE, ncol = nalt)
    probs_03[resp_test, task] <- exp(sum(log(tmp_XGamma) - log(sum(tmp_XGamma))))
  }
}

model_comparison[1, 4] <- round(mean(hits_01), 3)
model_comparison[1, 5] <- round(mean(probs_01), 3)
model_comparison[2, 4] <- round(mean(hits_02), 3)
model_comparison[2, 5] <- round(mean(probs_02), 3)
model_comparison[3, 4] <- round(mean(hits_03), 3)
model_comparison[3, 5] <- round(mean(probs_03), 3)

if (!file.exists(here::here("figures", str_c("model_comparison.rds")))) {
  full_model_comparison <- NULL
} else {
  full_model_comparison <- read_rds(here::here("figures", str_c("model_comparison.rds")))
}

full_model_comparison <- bind_rows(full_model_comparison, model_comparison)
write_rds(full_model_comparison, here::here("figures", str_c("model_comparison.rds")))

####################
# Compare to model_fit and predictive_fit_stacking.
####################

####################
# # Incorporate predictive_fit_ensemble.
# # predictive_fit_ensemble = function(indices, ensemble_weights, ensemble_draws, 
# #                                    test_X, test_Y, mat_ana, mat_screen, test_Z){
# #   # Compute the hit rate, hit prob, and loo metrics for the ensemble model.
# #   #   ensemble_weights - estimated weights for each of the models
# #   #   ensemble_draws - ensemble output with log_lik, betadraws, gammas, and Omegas for each model
# #   #   test_Y - choices (hold-out sample)
# #   #   test_X - design matrices (hold-out sample)
# #   #   test_Z - matrix of covariates
# #   
# #   nens <- length(ensemble_draws)
# #   ndraw <- length(ensemble_draws[[1]]$log_lik[,1,1]) # Number of draws
# #   nresp <- length(test_Y[,1])           # Number of respondents
# #   nscns <- length(test_X[1, ,1,1])      # Number of choice tasks
# #   nalts <- length(test_X[1,1, ,1])      # Number of alternatives 
# #   nlvls <- length(test_X[1,1,1, ])      # Number of att levels
# #   if( is.null(test_Z) ) test_Z <- matrix(1, nr = nresp, nc = 1)
# #   
# #   #weight log_lik for each model to get log_lik for ensemble
# #   LLmat_ens = matrix(0, nr=ndraw , 
# #                      # nc=exp(sum(log(dim(ensemble_draws[[k]]$log_lik))))/ndraw)
# #                      nc = dim(ensemble_draws[[1]]$log_lik)[2] * dim(ensemble_draws[[1]]$log_lik)[3])
# #   # loglik=ensemble_draws[[1]]$log_lik
# #   # ndraw=dim(loglik)[1]
# #   for(k in 1:nens){
# #     #extract log_lik array from each stanfit object
# #     loglik=ensemble_draws[[k]]$log_lik
# #     LLmat <- matrix(loglik, nr=ndraw)
# #     LLmat_ens <- LLmat_ens + ensemble_weights[k]*LLmat
# #   }
# #   
# #   #get effective sample size
# #   cores <- parallel::detectCores()
# #   r_eff_ens <- loo::relative_eff(x = exp(LLmat_ens), chain_id = double(ndraw)+1, cores = cores)
# #   
# #   #apply PSIS via loo to ensemble likelihoods (loo fit metrics)
# #   loo_fit_ens <- loo::loo.matrix(LLmat_ens, r_eff = r_eff_ens,
# #                                  cores = cores, save_psis = FALSE)
# #   
# #   # OR IS THERE A DEFAULT FOR THE EFFECTIVE SAMPLE SIZE?
# #   
# #   #stack resps and scns to avoid loops (this needs changed if using hold out tasks)
# #   test_X_stacked <- NULL
# #   for(resp in 1:nresp){
# #     for(scn in 1:nscns){
# #       test_X_stacked <- rbind(test_X_stacked,test_X[resp,scn,,])
# #     }
# #   }
# #   
# #   #stack scn choices to avoid loops
# #   test_Y_stacked <- matrix(t(test_Y),nc=1)
# #   
# #   #loop over ensemble models to calculate individual hit probs 
# #   probs_ens <- matrix(0, nr = nalts , nc = nresp*nscns)
# #   
# #   #loop over models
# #   for(model in 1:nens){
# #     #get betas for different hit rate calculations
# #     #using mean (over draws) of the mean of the post dist
# #     Umat <- matrix(0, nr = nresp*nscns*nalts)
# #     
# #     #get gammas
# #     meangammas=ensemble_draws[[model]]$Gamma
# #     
# #     # WHAT IS THIS?!
# #     # GAMMAS THAT AREN'T ATTENDED TO OR BEING SCREENED ON ARE UNIDENTIFIED
# #     # SO PASS ON THE INFORMATION FROM THE CLEVER RANDOMIZATION TO ZERO OUT
# #     # UNIDENTIFIED PARAMETERS FOR EACH ENSEMBLE MEMBER.
# #     
# #     # #adjust due to pathology approach in ensembles
# #     # index_ana <- mat_ana[model,]
# #     
# #     #loop over respondents
# #     for(resp in 1:nresp){
# #       #multiply by Z to get mean of dist of het for resp
# #       betas <- matrix(test_Z[resp,]%*%meangammas, nc=1)
# #       
# #       # WHAT IS THIS?!
# #       # #set gammadraws_mat column = 0 if ensemble ignores the level
# #       # betas[index_ana==1,] <- 0
# #       
# #       #get utility for each alternative
# #       Umat[((resp-1)*nalts*nscns+1):((resp)*nalts*nscns),] <- 
# #         exp(test_X_stacked[((resp-1)*nalts*nscns+1):((resp)*nalts*nscns),]%*%
# #               matrix(betas))
# #     }
# #     
# #     #find probabilities for each task and resp
# #     Umat_byscn <- matrix(Umat, nr = nalts) 
# #     sums <- t(matrix(rep(colSums(Umat_byscn),nalts), nc=nalts))
# #     #combine with other model probs weight by ensemble weights 
# #     probs_ens <- probs_ens + 
# #       (Umat_byscn/sums)*ensemble_weights[model]
# #   }
# #   
# #   #find location of highest prob meangammas
# #   locs <- apply(probs_ens,2,which.max)
# #   
# #   #calculate hits meangammas
# #   hits <- double(nresp*nscns)
# #   hits[locs==test_Y_stacked] <- 1
# #   
# #   #calculate hit probs
# #   hit_probs <- colSums(probs_ens*diag(nalts)[,test_Y_stacked])
# #   
# #   return(list(hit_rate = mean(hits), hit_prob = mean(hit_probs), loo_fit=loo_fit_ens))
# # }
# 
# predictive_fit_ensemble <- function(indices, ensemble_weights, ensemble_draws, test_X, test_Y, mat_ana, mat_screen, test_Z, ensemble_fit) {
#   # Compute the hit rate, hit prob, and loo metrics for the ensemble model.
#   #   ensemble_weights - estimated weights for each of the models
#   #   ensemble_draws - ensemble output with log_lik, betadraws, gammas, and Omegas for each model
#   #   test_Y - choices (hold-out sample)
#   #   test_X - design matrices (hold-out sample)
#   #   test_Z - matrix of covariates
#   
#   nmemb <- length(ensemble_draws) # Number of ensemble members.
#   nresp <- dim(test_X)[1]         # Number of respondents.
#   nscns <- dim(test_X)[2]         # Number of choice tasks.
#   
#   # Log-likelihood function for the MNL.
#   ll_mnl = function (beta, y, X) {
#     nvars = ncol(X)        # Number of attribute levels.
#     nscns = length(y)      # Number of choice tasks.
#     nalts = nrow(X)/nscns  # Number of alternatives.
#     
#     # Compute Xbeta across all choice tasks.
#     Xbeta = matrix(exp(X%*%beta),byrow=TRUE,ncol=nalts)
#     
#     # Numerator: Xbeta values associated with each choice.
#     choices = cbind(c(1:nscns),y)
#     numerator = Xbeta[choices]
#     
#     # Denominator: Xbeta values associated with each task.
#     iota = c(rep(1,nalts))
#     denominator = Xbeta%*%iota
#     
#     # Return the logit summed across choice tasks.
#     return(sum(log(numerator) - log(denominator)))
#   }
#   
#   # Compute LOO and predict Y.
#   memb_loo <- NULL
#   memb_hits <- NULL
#   memb_probs <- NULL
#   for (memb in 1:nmemb) {
#     # Get the mean of distribution of heterogeneity.
#     
#     ####################################################
#     # We use the average of Beta as the population mean without reference to Sigma.
#     ####################################################
#     
#     # Gamma_mean <- ensemble_draws[[memb]]$Gamma
#     # Gamma_mean <- ensemble_fit$ensemble_draws$Gamma |> 
#     #   select("mean") |> 
#     #   as.matrix()
#     Gamma_mean <- ensemble_draws[[memb]]$Beta |> 
#       spread_draws(Beta[, j]) |> 
#       summarize(mean = mean(Beta)) |> 
#       select("mean") |> 
#       as.matrix()
#     
#     # Predict Y.
#     hits <- NULL
#     probs <- NULL
#     for (resp in 1:nresp) {
#       for (scn in 1:nscns) {
#         # Pull the relevant Y and X.
#         Y_scn <- test_Y[resp, scn]
#         X_scn <- test_X[resp, scn,,]
#         
#         # Compute hit and probability.
#         # hits <- c(hits, Y_scn == which.max(X_scn %*% t(Gamma_mean)))
#         # probs <- c(probs, exp(ll_mnl(t(Gamma_mean), Y_scn, X_scn)))
#         
#         ####################################
#         # In-sample predictive fit using (modified) betas.
#         # Will need to save $Beta along with $Gamma.
#         ####################################
#         
#         hits <- c(hits, Y_scn == which.max(X_scn %*% Gamma_mean))
#         probs <- c(probs, exp(ll_mnl(Gamma_mean, Y_scn, X_scn)))
#       }
#     }
#     
#     # Compute LOO, the average hit rate, and hit probability.
#     # memb_loo <- c(memb_loo, loo(ensemble_fit[[memb]])$elpd_loo)
#     memb_loo <- c(memb_loo, 0)
#     memb_hits <- c(memb_hits, mean(hits))
#     memb_probs <- c(memb_probs, mean(probs))
#   }
#   
#   # Return the weighted average hit rate and hit probability.
#   return(list(
#     hit_rate = ensemble_weights %*% memb_hits, 
#     hit_prob = ensemble_weights %*% memb_probs, 
#     loo_fit = ensemble_weights %*% memb_loo
#   ))
# }
####################

####################
# # Compute fit metrics for ensemble.
# # ensemble_draws <- ensemble_fit$ensemble_draws$draws(format = "df", variables = c("Beta", "Gamma", "Omega", "tau"))
# # ensemble_pred_fit <- predictive_fit_hmnl(
# #   hmnl_draws = ensemble_fit$ensemble_draws, 
# #   test_X = data$test_X, 
# #   test_Y = data$test_Y,
# #   test_Z = data$test_Z
# # )
# ensemble_pred_fit <- predictive_fit_ensemble(
#   indices = c(ind_none, ind_ana, ind_screen, ind_ana_screen, ind_Z),
#   # ensemble_weights = ensemble_fit$ensemble_weights,
#   ensemble_weights = 1, # Equal weights.
#   ensemble_draws = ensemble_fit$ensemble_draws,
#   test_X = data$test_X,
#   test_Y = data$test_Y,
#   test_Z = data$test_Z,
#   mat_ana = ensemble_fit$mat_ana,
#   mat_screen = ensemble_fit$mat_screen,
#   ensemble_fit = ensemble_fit
# )
# 
# # Append results to the model comparison data frame.
# model_comparison <- model_comparison %>% 
#   bind_rows(
#     tibble(
#       Model = "Ensemble",
#       # LOO = ensemble_pred_fit$loo_fit$elpd_loo,
#       LOO = NA,
#       "Hit Rate" = ensemble_pred_fit$hit_rate,
#       "Hit Prob" = ensemble_pred_fit$hit_prob
#     )
#   )
# 
# # Print results.
# model_comparison
####################



########## # Temp commented.
# # Try equal weights.
# ensemble_fit$ensemble_weights <- rep(NA, length(ensemble_fit$ensemble_draws))
# for (k in 1:length(ensemble_fit$ensemble_weights)) {
#   ensemble_fit$ensemble_weights[k] <- 1 / length(ensemble_fit$ensemble_weights)
# }
# 
# # # Try random weights.
# # temp_ensemble_weights <- runif(n = length(ensemble_fit$ensemble_draws), min = 0, max = 1)
# # for (k in 1:length(ensemble_fit$ensemble_weights)) {
# #   ensemble_fit$ensemble_weights[k] <- temp_ensemble_weights[k] / sum(temp_ensemble_weights)
# # }
# 
# # # Use the first half of the test data for validation data.
# # data$validate_Y <- data$test_Y[1:round(dim(data$test_Y)[1]/2),]
# # data$validate_X <- data$test_X[1:round(dim(data$test_X)[1]/2),,,]
# # data$validate_Z <- as.matrix(data$test_Z[1:round(dim(data$test_Z)[1]/2),])
# # data$test_Y <- data$test_Y[(round(dim(data$test_Y)[1]/2) + 1):dim(data$test_Y)[1],]
# # data$test_X <- data$test_X[(round(dim(data$test_X)[1]/2) + 1):dim(data$test_X)[1],,,]
# # data$test_Z <- as.matrix(data$test_Z[(round(dim(data$test_Z)[1]/2) + 1):dim(data$test_Z)[1],])
# # 
# # # Produce predictions for each ensemble member using the validation data.
# # ensemble_predictions <- vector(mode = "list", length = nmember)
# # for (k in 1:nmember) {
# #   ensemble_predictions[[k]] <- predictive_fit_stacking(
# #     member_draws = ensemble_fit$ensemble_draws[[k]],
# #     validate_X = data$validate_X,
# #     validate_Z = data$validate_Z
# #   )
# # }
# # 
# # # Produce counts or probabilities to calculate the ensemble weights.
# # meta_Y <- as.vector(t(data$validate_Y))
# # meta_pred_X <- matrix(NA, nrow = length(meta_Y), ncol = nmember)
# # meta_prob_X <- array(NA, dim = c(length(meta_Y), max(meta_Y), nmember))
# # for (k in 1:nmember) {
# #   for (n in 1:length(meta_Y)) {
# #     meta_pred_X[n,k] <- ensemble_predictions[[k]]$predicted_Y[n]
# #     meta_prob_X[n,,k] <- ensemble_predictions[[k]]$predicted_probs[,n]
# #   }
# # }
# # temp_ensemble_counts <- rep(NA, nmember)
# # temp_ensemble_sum_probs <- rep(NA, nmember)
# # for(k in 1:nmember) {
# #   temp_ensemble_counts[k] <- sum(meta_Y == meta_pred_X[,k])
# #   temp_probs <- NULL
# #   for (n in 1:length(meta_Y)) {
# #     temp_probs <- c(temp_probs, meta_prob_X[n, meta_Y[n], k])
# #   }
# #   temp_ensemble_sum_probs[k] <- sum(temp_probs)
# # }
# # 
# # # Normalize the counts or probabilities.
# # for (k in 1:nmember) {
# #   ensemble_fit$ensemble_weights[k] <- temp_ensemble_counts[k] / sum(temp_ensemble_counts)
# #   # ensemble_fit$ensemble_weights[k] <- temp_ensemble_sum_probs[k] / sum(temp_ensemble_sum_probs)
# # }
# 
# # # Restructure validation data and predictions or probabilities for the meta-learner.
# # meta_Y <- as.vector(t(data$validate_Y))
# # meta_X <- array(NA, dim = c(length(meta_Y), max(meta_Y), nmember))
# # # for (n in 1:length(meta_Y)) {
# # #   temp_X <- NULL
# # #   for (k in 1:nmember) {
# # #     temp_X <- cbind(temp_X, as.vector(t(ensemble_predictions[[k]]$predicted_Y))[n])
# # #   }
# # #   meta_X[n,,] <- matrix(temp_X, nrow = max(meta_Y), byrow = TRUE)
# # # }
# # for (k in 1:nmember) {
# #   for (n in 1:length(meta_Y)) {
# #     meta_X[n,,k] <- ensemble_predictions[[k]]$predicted_probs[,n]
# #   }
# # }
# # 
# # # Produce weights for each of the choice tasks in the validation data.
# # stan_data <- list(
# #   N = dim(meta_X)[1], # Number of observations.
# #   A = dim(meta_X)[2], # Number of choice alternatives.
# #   L = dim(meta_X)[3], # Number of (estimable) attribute levels.
# # 
# #   Y = meta_Y,         # Vector of observations.
# #   X = meta_X          # Matrix of observation-level covariates.
# # )
# # 
# # options(mc.cores = parallel::detectCores())
# # rstan_options(auto_write = FALSE)
# # 
# # meta_fit <- stan(
# #   here::here("Code", "src", "mnl.stan"),
# #   data = stan_data,
# #   seed = 42
# # )
# # 
# # # Extract weights for each ensemble and normalize.
# # temp_ensemble_weights <- extract(meta_fit, pars = c("beta"))
# # temp_ensemble_weights <- apply(temp_ensemble_weights$beta, 2, mean)
# # temp_ensemble_weights <- (temp_ensemble_weights - min(temp_ensemble_weights)) /
# #   (max(temp_ensemble_weights) - min(temp_ensemble_weights))
# # for (k in 1:nmember) {
# #   ensemble_fit$ensemble_weights[k] <- temp_ensemble_weights[k] / sum(temp_ensemble_weights)
# # }
# 
# # Check the ensemble weights.
# sum(ensemble_fit$ensemble_weights)
# hist(ensemble_fit$ensemble_weights)
# 
# ensemble_fit$ensemble_weights
# 
# # # Try dropping the final member and normalizing weights.
# # ensemble_fit$ensemble_weights[length(ensemble_fit$ensemble_weights)] <- 0
# # for (k in 1:length(ensemble_fit$ensemble_weights)) {
# #   ensemble_fit$ensemble_weights[k] <- (ensemble_fit$ensemble_weights[k] / sum(ensemble_fit$ensemble_weights))
# # }
# 
# # # Try weeding out ensemble members with VB errors.
# # ensemble_fit$ensemble_draws <-
# #   ensemble_fit$ensemble_draws[-which(lapply(ensemble_fit$ensemble_draws, length) == 1)]
# # length(ensemble_fit$ensemble_draws)
# 
# # Compute Model Fit -------------------------------------------------------
# # Extract needed draws.
# hmnl_draws <- extract(hmnl_fit, pars = c("Beta", "Gamma", "Omega", "tau"))
# 
# # Compute HMNL predictive fit.
# hmnl_pred_fit <- predictive_fit_hmnl(
#   hmnl_draws = hmnl_draws, 
#   test_X = data$test_X, 
#   test_Y = data$test_Y,
#   test_Z = data$test_Z
# )
# 
# # Create a model comparison data frame.
# model_comparison <- tibble(
#   Model = "HMNL",
#   LOO = loo(hmnl_fit)$elpd_loo,
#   "Hit Rate" = hmnl_pred_fit$hit_rate,
#   "Hit Prob" = hmnl_pred_fit$hit_prob
# )
# 
# # Print results.
# model_comparison
# 
# # Compute fit metrics for ensemble.
# ensemble_pred_fit <- predictive_fit_ensemble(
#   indices = c(ind_none, ind_ana, ind_screen, ind_ana_screen, ind_Z),
#   ensemble_weights = ensemble_fit$ensemble_weights, 
#   ensemble_draws = ensemble_fit$ensemble_draws, 
#   test_X = data$test_X, 
#   test_Y = data$test_Y,
#   test_Z = data$test_Z,
#   mat_ana = ensemble_fit$mat_ana,
#   mat_screen = ensemble_fit$mat_screen,
#   ensemble_fit = ensemble_fit
# )
# 
# # Append results to the model comparison data frame.
# model_comparison <- model_comparison %>% 
#   bind_rows(
#     tibble(
#       Model = "Ensemble",
#       # LOO = ensemble_pred_fit$loo_fit$elpd_loo,
#       LOO = NA,
#       "Hit Rate" = ensemble_pred_fit$hit_rate,
#       "Hit Prob" = ensemble_pred_fit$hit_prob
#     )
#   )
# 
# # Print results.
# model_comparison
# 
# # Still need to add competing models here with predictive fit.
# 
# # Compute fit metrics for ANA model
# #ana_pred_fit <- predictive_fit_ana(
# #  ana_fit = ana_draws, 
# #  test_X = data$test_X, 
# #  test_Y = data$test_Y,
# #  Z=NULL
# #)
# 
# # Append results to the model comparison data frame.
# model_comparison <- model_comparison %>% 
#   bind_rows(
#     tibble(
#       Model = "ANA",
#       LOO = NA,
#       "Hit Rate" = NA,
#       "Hit Prob" = NA
#     )
#   )
# 
# #ConjScreen
# # Compute fit metrics for Conj Screening model
# #cscreen_pred_fit <- predictive_fit_cscreen(
# #  cscreen_fit = cscreen_draws, 
# #  test_X = data$test_X, 
# #  test_Y = data$test_Y,
# #  Z=NULL
# #)
# 
# # Append results to the model comparison data frame.
# model_comparison <- model_comparison %>% 
#   bind_rows(
#     tibble(
#       Model = "Conj Screen",
#       LOO = NA,
#       "Hit Rate" = NA,
#       "Hit Prob" = NA
#     )
#   )
# 
# # # Save model comparison data frame.
# # write_rds(model_comparison, here::here("Figures", str_c("model-fit_", file_id, "_", nmember, ".rds")))
# 
# # Parameter Recovery ------------------------------------------------------
# if (ind_test == 1) {
#   # Compare parameter estimates and the true values.
#   tibble(
#     param = as.factor(1:length(data$bbar)),
#     true = data$bbar,
#     hmnl_estimate = apply(hmnl_draws$Gamma, 3, mean),
#     ensm_estimate_1 = as.vector(ensemble_fit$ensemble_draws[[1]]$Gamma),
#     ensm_estimate_2 = as.vector(ensemble_fit$ensemble_draws[[2]]$Gamma)
#   ) %>% 
#     mutate(
#       hmnl_est_diff = true - hmnl_estimate,
#       ensm_est_diff_1 = true - ensm_estimate_1,
#       ensm_est_diff_2 = true - ensm_estimate_2
#     ) %>% 
#     pivot_longer(contains("diff"), names_to = "contrast", values_to = "difference") %>%
#     ggplot(aes(x = param, y = difference, fill = contrast)) +
#     geom_col(position = "dodge") +
#     scale_fill_discrete(type = c("gray40", "gray50", "darkred")) +
#     # ylim(-2, 2) +
#     theme(legend.position = "bottom")
#   
#   ggsave(here::here("Figures", str_c("parameter-recovery_", file_id, ".png")), width = 7, height = 3)
# }
# 
