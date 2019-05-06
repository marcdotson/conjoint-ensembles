// Hierarchical multinomial logit with a multivariate normal distribution of 
// heterogeneity for modeling choice data from a conjoint experiment.
data {
  int<lower=1> N;                 // Number of respondents.
  int<lower=1> S;                 // Number of choice tasks per respondent.
  int<lower=2> P;                 // Number of product alternatives per choice task.
  int<lower=1> L;                 // Number of (estimable) attribute levels.
  int<lower=1> C;                 // Number of respondent-level covariates.
  
  int<lower=1, upper=P> Y[N, S];  // Matrix of observed choices.
  matrix[P, L] X[N, S];           // Array of experimental designs per choice task.
  matrix[N, C] Z;                 // Matrix of respondent-level covariates.
  
  real mu_loc;                    // Location of the means of Beta.
  real<lower=0> mu_scale;         // Scale of the means of Beta.
  real alpha_loc;                 // Location of the variance of Beta.
  real<lower=0> alpha_scale;      // Scale of the variance of Beta.
  real<lower=0> lkj_corr_shape;   // Correlation matrix hyperprior.
}

parameters {
  matrix[L, N] alpha;                        // Prior on variance of Beta.
  cholesky_factor_corr[L] L_Omega;           // Cholesky factorization of hyperprior.
  vector<lower=0, upper=pi()/2>[L] tau_unif; // ??
  matrix[C, L] mu;                           // Prior on mean of Beta.
}

transformed parameters {
  matrix[N, L] Beta;                                         // Matrix of Beta coefficients.
  vector<lower=0>[L] tau;                                    // Prior scale.
  for (l in 1:L) tau[l] = 2.5 * tan(tau_unif[l]);            // ??
  Beta = Z * mu + (diag_pre_multiply(tau,L_Omega) * alpha)'; // ??
}

model {
  //priors
  to_vector(alpha) ~ normal(alpha_loc, alpha_scale);
  L_Omega ~ lkj_corr_cholesky(lkj_corr_shape);
  to_vector(mu) ~ normal(mu_loc, mu_scale);

  // model fitting
  for (r in 1:N) {
    for (t in 1:S) {
      Y[r,t] ~ categorical_logit(X[r,t]*Beta[r]');
    }
  }
}

generated quantities {
  // Yp is predicted choices for new data.
  real Y_ppc[N, S];
  vector[N*S] log_lik;
  {
    matrix[N, S] temp_log_lik;
    for (r in 1:N) {
      for (t in 1:S) {
        Y_ppc[r, t] = categorical_logit_rng(X[r, t] * Beta[r]');
        temp_log_lik[r, t] = categorical_logit_lpmf(Y[r, t] | X[r, t] * Beta[r]');
      }
    }
    log_lik = to_vector(temp_log_lik);
  }
}
