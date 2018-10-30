// HBMNL for ANA pathology
data {
  int<lower=2> A; // number of alternatives (choices) per question
  int<lower=1> L; // number of feature variables
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
  real<lower=0> ksi; // chi-squared for multi student t
}

transformed parameters {
  matrix[R, L] B; // matrix of beta coefficients
  vector<lower=0>[L] tau; // prior scale
  real<lower=0> nu = 1.0; // prior for multi student t distribution

  for (l in 1:L) tau[l] = 2.5 * tan(tau_unif[l]);
  B = Z * mu + nu*ksi*(diag_pre_multiply(tau,L_Omega) * alpha)';
}

model {
  //priors
  to_vector(alpha) ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(.5);
  to_vector(mu) ~ normal(0, 1);
  ksi ~ inv_chi_square(nu);

  // model fitting
  for (r in 1:R) {
    for (t in 1:T) {
      Y[r,t] ~ categorical_logit(X[r,t]*B[r]');
    }
  }
}

generated quantities {
  // Yp is predicted choices for new data.
  real Y_ppc[R, T];
  for (r in 1:R) {
    for (t in 1:T) {
      Y_ppc[r,t] = categorical_logit_rng(X[r,t] * B[r]');
    }
  }
}
