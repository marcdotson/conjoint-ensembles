// HBMNL for discrete choice experiments
data {
  int<lower=2> A; // number of alternatives (choices) per question
  int<lower=1> L; // number of feature variables
  int<lower=1> R; // number of respondents
  int<lower=1> T; // number of questions (unique inquiries)
  int<lower=1> C; // number of respondent covariates (demographics, etc)
  int<lower=1> Rtest;
  int<lower=1> Ttest;

  matrix[A, L] X[R, T]; // matrix of attributes for each obs
  int<lower=1, upper=A> Y[R, T]; // observed responses
  matrix[R, C] Z; // vector of covariates for each respondent

  matrix[A, L] Xtest[Rtest, Ttest]; // test design matrix

  real mu_mean;
  real alpha_mean;
  real<lower=0> mu_scale;
  real<lower=0> alpha_scale;
  real<lower=0> lkj_param;
}

parameters {
  matrix[L, R] alpha; // prior on variance of utilities B
  vector<lower=0, upper=pi()/2>[L] tau_unif;
  matrix[C, L] mu; // prior on mean of utilities B
  cholesky_factor_corr[L] L_Omega;

}

transformed parameters {
  matrix[R, L] B; // matrix of beta coefficients
  vector<lower=0>[L] tau; // prior scale
  for (l in 1:L) tau[l] = 2.5 * tan(tau_unif[l]);
  B = Z * mu + (diag_pre_multiply(tau,L_Omega) * alpha)';
}

model {
  //priors
  to_vector(alpha) ~ normal(alpha_mean, alpha_scale);
  to_vector(mu) ~ uniform(mu_mean, mu_scale);
  L_Omega ~ lkj_corr_cholesky(lkj_param);


  // model fitting
  for (r in 1:R) {
    for (t in 1:T) {
      Y[r, t] ~ categorical_logit(X[r, t]*B[r]');
    }
  }
}

generated quantities {
  // Yp is predicted choices for new data.
  int<lower=0, upper=A> Yhat[Rtest, Ttest];
  int<lower=0, upper=1> Yc[Rtest, Ttest, A];
  vector[L] B_avg;

  for (r in 1:Rtest) {
    for (t in 1:Ttest) {
      for (l in 1:L) B_avg[l] = mean(B[,l]);
      Yhat[r, t] = categorical_logit_rng(Xtest[r, t]*B_avg);
      for (a in 1:A) {
        if (Yhat[r, t] == a) {
          Yc[r, t, a] = 1;
        }
        else {
          Yc[r, t, a] = 0;
        }
      }
    }
  }
}
