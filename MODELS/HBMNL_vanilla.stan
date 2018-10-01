// HBMNL for discrete choice experiments
data {
  int<lower=2> C; // number of alternatives (choices) per question
  int<lower=1> K; // number of feature variables
  int<lower=1> J; // number of respondents
  int<lower=1> S; // number of questions (unique inquiries)
  int<lower=1> G; // number of respondent covariates (demographics, etc)
  int<lower=1,upper=C> Y[J, S]; // observed responses
  matrix[C, K] X[J, S]; // matrix of attributes for each obs
  matrix[J, G] Z; // vector of covariates for each respondent
}

parameters {
  matrix[K, J] alpha; // prior on variance of utilities B
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0,upper=pi()/2>[K] tau_unif;
  matrix[G, K] mu; // prior on mean of utilities B
  real<lower=0> sigma;
}

transformed parameters {
  matrix[J, K] B; // matrix of beta coefficients
  vector<lower=0>[K] tau; // prior scale
  for (k in 1:K) tau[k] = 2.5 * tan(tau_unif[k]);
  B = Z * mu + (diag_pre_multiply(tau,L_Omega) * alpha)';
}

model {
  //priors
  to_vector(alpha) ~ normal(0, 10);
  L_Omega ~ lkj_corr_cholesky(5);
  to_vector(mu) ~ normal(0, 1);

  // model fitting
  for (j in 1:J) {
    for (s in 1:S) {
      Y[j,s] ~ categorical_logit(X[j,s]*B[j]');
    }
  }
}

generated quantities {
  // Yp is predicted choices for new data.
  real log_lik[J, S];
  for (j in 1:J) {
    for (s in 1:S) {
      log_lik[j,s] = categorical_logit_lpmf( Y[j,s] | to_vector(B[j]) );
    }
  }
}
