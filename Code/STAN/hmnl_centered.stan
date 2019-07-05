// Data values, hyperparameters, observed choices, and the experimental design.
data {
  int<lower=1> R; // Number of respondents.
  int<lower=1> T; // Number of choice tasks per respondent.
  int<lower=2> A; // Number of product alternatives per choice task.
  int<lower=1> L; // Number of (estimable) attribute levels.
  int<lower=1> C; // Number of respondent-level covariates.
  int<lower=1> Rtest;
  int<lower=1> Ttest;
  

  matrix[A, L] X[R, T];          // Array of experimental designs per choice task.
  int<lower=1, upper=A> Y[R, T]; // Matrix of observed choices.
  matrix[R, C] Z;                // Matrix of respondent-level covariates.

  matrix[A, L] Xtest[Rtest, Ttest]; // test design matrix

  real mu_mean;           // Mean of coefficients for the heterogeneity model.
  real<lower=0> mu_scale; // Scale of coefficients for the heterogeneity model.
  real tau_mean;             // Mean of scale parameters for the heterogeneity model.
  real<lower=0> tau_scale;   // Scale of scale parameters for the heterogeneity model.
  real<lower=0> Omega_shape; // Shape of correlation matrix for the heterogeneity model.
}

parameters {
  matrix[R, L] B;        // Matrix of beta (part-worth) coefficients.
  matrix[C, L] mu;       // Matrix of coefficients for the heterogeneity model.
  vector<lower = 0>[L] tau; // Vector of scale parameters for the heterogeneity model.
  corr_matrix[L] Omega;     // Correlation matrix for the heterogeneity model.
}

transformed parameters {
  // Covariance matrix for the heterogeneity model.
  cov_matrix[L] Sigma = quad_form_diag(Omega, tau);
  matrix[R, L] Zmu = Z * mu;
}

model {
  // Hyperpriors on mu, tau, and Omega (and thus Sigma).
  to_vector(mu) ~ normal(mu_mean, mu_scale);
  tau ~ cauchy(tau_mean, tau_scale);
  Omega ~ lkj_corr(Omega_shape);
  
  // Hierarchical multinomial logit.
  for (r in 1:R) {
    B[r, ] ~ multi_normal(Z[r,] * mu, Sigma);    
    for (t in 1:T) {
      Y[r, t] ~ categorical_logit(X[r, t] * B[r,]');
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
