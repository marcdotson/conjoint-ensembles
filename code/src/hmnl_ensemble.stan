// Index values, hyperprior values, observations, and covariates.
data {
  int<lower = 1> R;               // Number of respondents.
  int<lower = 1> S;               // Number of choice tasks.
  int<lower = 2> A;               // Number of choice alternatives.
  int<lower = 1> I;               // Number of observation-level covariates.
  int<lower = 1> J;               // Number of population-level covariates.

  real Gamma_mean;                // Mean of population-level means.
  real<lower=0> Gamma_scale;      // Scale of population-level means.
  real<lower=0> Omega_shape;      // Shape of population-level scale.
  real tau_mean;                  // Mean of population-level scale.
  real<lower=0> tau_scale;        // Scale of population-level scale.

  array[R, S] int Y;              // Array of observations.
  array[R, S] matrix[A, I] X;     // Array of observation-level covariates.
  matrix[R, J] Z;                 // Matrix of population-level covariates.
  array[R, I] int array_ana;      // Array of ensemble indicators for ANA.
  array[R, I] int array_screen;   // Array of ensemble indicators for screening.
  array[R, 1] int array_qual;     // Array of ensemble indicators for respondent quality.
}

// Parameters and hyperparameters.
parameters {
  matrix[J, I] Gamma;                // Matrix of population-level hyperparameters.
  corr_matrix[I] Omega;              // Population model correlation matrix hyperparameters.
  vector<lower = 0>[I] tau;          // Population model vector of scale hyperparameters.
  matrix[R, I] Delta;                // Matrix of non-centered observation-level parameters.
}

// Deterministic transformation.
transformed parameters {
  // Matrix of centered observation-level parameters.
  matrix[R, I] Beta;
  
  // Impose fixed values using ANA indicator array.
  for (r in 1:R) {
    for (i in 1:I) {
      if (array_ana[r, i] == 1) {
        Beta[r, i] = 0;
      } else {
        Beta[r, i] = Beta[r, i];
      }
    }
  }
  
  // Impose fixed values using screening indicator array.
  for (r in 1:R) {
    for (i in 1:I) {
      if (array_screen[r, i] == 1) {
        Beta[r, i] = -100;
      } else {
        Beta[r, i] = Beta[r, i];
      }
    }
  }
  
  // Impose fixed values using respondent quality indicator array.
  for (r in 1:R) {
    if (array_qual[r, 1] == 1) {
      Beta[r,] = Beta[r,] * 0;
    } else {
      Beta[r,] = Beta[r,];
    }
  }
  
  // Non-centered parameterization.
  for (r in 1:R) {
    Beta[r,] = Z[r,] * Gamma + Delta[r,] * quad_form_diag(Omega, tau);
  }
}

// Hierarchical multinomial logit model.
model {
  // Hyperpriors.
  to_vector(Gamma) ~ normal(Gamma_mean, Gamma_scale);
  Omega ~ lkj_corr(Omega_shape);
  tau ~ normal(tau_mean, tau_scale);

  // Non-centered population model and likelihood.
  for (r in 1:R) {
    Delta[r,] ~ normal(0, 1);
    for (s in 1:S) {
      Y[r, s] ~ categorical_logit(X[r, s,,] * Beta[r,]');
    }
  }
}

// Generated quantities conditioned on parameter draws.
generated quantities {
  matrix[R, S] log_lik; // Matrix of log likelihood values.
  matrix[I, I] Sigma;   // Covariance matrix for the population model.
  for (r in 1:R) {
    for (s in 1:S) {
      log_lik[r, s] = categorical_logit_lpmf(Y[r, s] | X[r, s,,] * Beta[r,]');
    }
  }
  Sigma = quad_form_diag(Omega, tau);
}
