// HBMNL for discrete choice experiments
data {
  int<lower=2> A; // number of alternatives (choices) per question
  int<lower=1> L; // number of feature variables (levels)
  int<lower=1> R; // number of respondents
  int<lower=1> T; // number of questions (unique inquiries)
  int<lower=1> C; // number of respondent covariates (demographics, etc)
  int<lower=1,upper=A> Y[R, T]; // observed responses
  matrix[A, L] X[R, T]; // matrix of attributes for each obs
  matrix[R, C] Z; // vector of covariates for each respondent
}

parameters {
  matrix[L, R] alpha; // prior on variance of utilities B
  cholesky_factor_corr[L] L_Omega;
  vector<lower=0,upper=pi()/2>[L] tau_unif;
  matrix[C, L] mu; // prior on mean of utilities B
  vector[L] ksi; // bernoulli process prior representing ANA pathology
  real<lower=0, upper=1> theta[L]; // probabilities associated with bernoulli process
}

transformed parameters {
  matrix[R, L] B; // matrix of beta coefficients
  vector<lower=0>[L] tau; // prior scale
  for (k in 1:L) tau[k] = 2.5 * tan(tau_unif[k]);
  B = Z * mu + (diag_pre_multiply(tau,L_Omega) * alpha)';
}

model {
  //priors
  to_vector(alpha) ~ normal(0, 10);
  L_Omega ~ lkj_corr_cholesky(5);
  to_vector(mu) ~ normal(0, 1);
  ksi ~ bernoulli(theta);

  // model fitting
  for (r in 1:R) {
    for (t in 1:T) {
      Y[r,t] ~ categorical_logit(X[r,t]*B[r]');
    }
  }
}

generated quantities {
  // Yp is predicted choices for new data.
  real log_lik[R, T];
  for (r in 1:R) {
    for (t in 1:T) {
      log_lik[r,t] = categorical_logit_lpmf( Y[r,t] | to_vector(B[r]) );
    }
  }
}
