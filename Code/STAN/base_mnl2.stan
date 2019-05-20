// HBMNL for discrete choice experiments
data {
  int<lower=2> A; // number of alternatives (choices) per question
  int<lower=1> L; // number of feature variables
  int<lower=1> R; // number of respondents
  int<lower=1> T; // number of questions (unique inquiries)
  int<lower=1> Rtest;
  int<lower=1> Ttest;

  matrix[A, L] X[R, T]; // matrix of attributes for each obs
  int<lower=1, upper=A> Y[R, T]; // observed responses

  matrix[A, L] Xtest[Rtest, Ttest]; // test design matrix

}

transformed data {
  int<lower=1> C = 1;
  matrix[R, C] Z = rep_matrix(1, R, C); // vector of covariates for each respondent
  real mu_mean = 0;
  real mu_scale = 1;
  real<lower=0> alpha_mean = 0;
  real<lower=0> alpha_scale = 10;
  real<lower=0> lkj_param = 5;
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
  to_vector(mu) ~ normal(mu_mean, mu_scale);
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
