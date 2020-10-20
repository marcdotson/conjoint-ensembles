// Index values, hyperprior values, observations, and covariates.
data {
  int<lower = 1> R;                  // Number of respondents.
  int<lower = 1> S;                  // Number of choice tasks.
  int<lower = 2> A;                  // Number of choice alternatives.
  int<lower = 1> I;                  // Number of observation-level covariates.
  int<lower = 1> J;                  // Number of population-level covariates.

  real Gamma_mean;                   // Mean of population-level means.
  real<lower=0> Gamma_scale;         // Scale of population-level means.
  real<lower=0> Omega_shape;         // Shape of population-level scale.
  real tau_mean;                     // Mean of population-level scale.
  real<lower=0> tau_scale;           // Scale of population-level scale.

  int<lower = 1, upper = A> Y[R, S]; // Matrix of observations.
  matrix[A, I] X[R, S];              // Array of observation-level covariates.
  matrix[R, J] Z;                    // Matrix of population-level covariates.
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

  // Non-centered parameterization.
  for (r in 1:R) {
    Beta[r,] = Z[r,] * Gamma + Delta[r,] * quad_form_diag(Omega, tau);
  }
  
  // Imposing clever randomization.
  // Either loop through all elements to find which are set to 0 (?) or determine
  // this beforehand and input as part of the data block, then reference those 
  // specific elements and set them to 0 (?) with something like:
  // Beta[, i] = 0;
  // 
  // Something like
  // real nondiag_xi[N, S, S];
  // nondiag_xi = xi;
  // for (s in 1:S)
  //   for (n in 1:N)
  //     nondiag_xi[n, S, S] = 0.0;
  // 
  // Or:
  // matrix[4, 6] x1[2, 1, 7];
  // matrix[4, 6] x2[3, 1, 7];
  // matrix[4, 6] x3[5, 1, 7];
  // x3 = append_array(x1, x2);
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
      Y[r, s] ~ categorical_logit(X[r, s] * Beta[r,]');
    }
  }
}

// Generated quantities conditioned on parameter draws.
generated quantities {
  // Compute log likelihood for model fit.
  matrix[R, S] log_lik;
  for (r in 1:R) {
    for (s in 1:S) {
      log_lik[r, s] = categorical_logit_lpmf(Y[r, s] | X[r, s] * Beta[r,]');
    }
  }
}
