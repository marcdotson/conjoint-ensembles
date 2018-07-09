data {
  int<lower = 2> nresp;                                 // Number of respondents.
  int<lower = 1> nscns;                                 // Number of choice scenarios for each respondent.
  int<lower = 1> nalts;                                 // Number of alternatives in each choice scenario.
  int<lower = 1> nlvls;                                 // Number of attribute levels for each alternative.
  int<lower = 1> ncovs;                                 // Number of covariates for each respondent.
  int<lower = 1, upper = nalts> Y[nresp, nscns];        // Picked alternatives for each respondent and choice scenario.
  matrix[nalts, nlvls] X[nresp, nscns];                 // Array of design matrices for each respondent and choice scenario.
  matrix[ncovs, nresp] Z;                               // Matrix of covariates for each respondent.
}
parameters {
  matrix[nlvls, nresp] Beta;                            // Respondent-level part-worths.
  matrix[nlvls, ncovs] Gamma;                           // Aggregate-level mean preference heterogeneity.
  corr_matrix[nlvls] Omega;                             // Aggregate-level correlation matrix of preference heterogeneity.
  vector<lower=0>[nlvls] tau;                           // Hyperparameter for reparameterization of Omega.
}
transformed parameters {
  cov_matrix[nlvls] Vbeta = quad_form_diag(Omega, tau); // Reparameterization of the Omega and tau to produce Vbeta.
}
model {
  // Priors.
  to_vector(Gamma) ~ normal(0, 10);
  tau ~ cauchy(0, 2.5);
  Omega ~ lkj_corr(2);
  
  // Likelihood.
  for (resp in 1:nresp) {
    // Draw from the distribution of heterogeneity.	
    Beta[, resp] ~ multi_normal(Gamma * Z[, resp], Vbeta);
    // Multionomial logit.
    for (scn in 1:nscns) {
      Y[nresp, nscns] ~ categorical_logit(X[resp, scn] * Beta[, resp]);
    }
  }
}
