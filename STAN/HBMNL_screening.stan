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
  vector<lower=0>[A] P[R];
}

transformed parameters {
  matrix[R, L] B; // matrix of beta coefficients
  vector<lower=0>[L] tau; // prior scale

  for (l in 1:L) tau[l] = 2.5 * tan(tau_unif[l]);
  B = Z * mu + (diag_pre_multiply(tau,L_Omega) * alpha)';
}

model {
  //priors
  to_vector(alpha) ~ normal(0, 10);
  L_Omega ~ lkj_corr_cholesky(2);
  to_vector(mu) ~ normal(0, 1);

  // model fitting
  for (r in 1:R) {
    P[r] ~ inv_gamma(0.5, 0.5);
    for (t in 1:T) {
      Y[r, t] ~ categorical_logit(X[r,t]*B[r]' - P[r]);
    }
  }
}

generated quantities {
  // Yp is predicted choices for new data.
  real Y_ppc[R, T];
  for (r in 1:R) {
    for (t in 1:T) {
      Y_ppc[r,t] = categorical_logit_rng(X[r,t] * B[r]' - P[r]);
    }
  }
}
