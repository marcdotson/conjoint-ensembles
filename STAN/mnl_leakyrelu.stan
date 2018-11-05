// HBMNL for discrete choice experiments
functions {
  matrix LeakyReLU(matrix X) {
    matrix[rows(X), cols(X)] X_new;
    for (j in 1:cols(X)) {
      for (i in 1:rows(X)) {
        if (X[i, j] <= -1) X_new[i, j] = -square(X[i, j]);
        else X_new[i, j] = X[i, j];
      }
    }
    return X;
  }
}

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
}

transformed parameters {
  matrix[R, L] B; // matrix of beta coefficients
  vector<lower=0>[L] tau; // prior scale
  for (l in 1:L) tau[l] = 2.5 * tan(tau_unif[l]);
  B = Z * mu + (diag_pre_multiply(tau,L_Omega) * alpha)';
  B = LeakyReLU(B);
}

model {
  //priors
  to_vector(alpha) ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(2);
  to_vector(mu) ~ normal(0, 5);

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
  matrix[R, T] log_lik;
  for (r in 1:R) {
    for (t in 1:T) {
      Y_ppc[r, t] = categorical_logit_rng(X[r, t] * B[r]');
      log_lik[r, t] = categorical_logit_lpmf(Y[r, t] | X[r, t] * B[r]');
    }
  }
}
