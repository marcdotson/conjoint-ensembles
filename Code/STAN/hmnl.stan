// hierarchical multinomial logit for discrete choice experiments
data {
  int<lower=2> A; // number of alternatives (choices) per question
  int<lower=1> L; // number of feature variables
  int<lower=1> R; // number of respondents
  int<lower=1> T; // number of questions (unique inquiries)
  int<lower=1> C; // number of respondent covariates (demographics, etc)
  int<lower=1, upper=A> Y[R, T]; // observed responses
  matrix[A, L] X[R, T]; // matrix of attributes for each obs
  matrix[R, C] Z; // vector of covariates for each respondent
  real mu_loc; // location of the means of B
  real<lower=0> mu_scale; // scale of the means of B
  real alpha_loc; // location of the variance of B
  real<lower=0> alpha_scale; // scale of the variance of B
  real<lower=0> lkj_corr_shape; // for correlation matrix hyperprior
}

parameters {
  matrix[L, R] alpha; // prior on variance of utilities B
  cholesky_factor_corr[L] L_Omega;
  vector<lower=0, upper=pi()/2>[L] tau_unif;
  matrix[C, L] mu; // prior on mean of utilities B
}

transformed parameters {
  matrix[R, L] B; // matrix of beta coefficients
  vector<lower=0>[L] tau; // prior scale
  for (l in 1:L) tau[l] = 2.5 * tan(tau_unif[l]);
  B = Z * mu + (diag_pre_multiply(tau,L_Omega) * alpha)';
}

model {
  //priors
  to_vector(alpha) ~ normal(alpha_loc, alpha_scale);
  L_Omega ~ lkj_corr_cholesky(lkj_corr_shape);
  to_vector(mu) ~ normal(mu_loc, mu_scale);

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
  vector[R*T] log_lik;
  {
    matrix[R, T] temp_log_lik;
    for (r in 1:R) {
      for (t in 1:T) {
        Y_ppc[r, t] = categorical_logit_rng(X[r, t] * B[r]');
        temp_log_lik[r, t] = categorical_logit_lpmf(Y[r, t] | X[r, t] * B[r]');
      }
    }
    log_lik = to_vector(temp_log_lik);
  }
}
