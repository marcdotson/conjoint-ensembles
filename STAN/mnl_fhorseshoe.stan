// Finnish Horseshoe prior for ANA
// Adapted from Michael Betancourt's case study on Bayes Sparse Regression
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

transformed data {
  real m0 = 0.5*L;
  real sig = 2; // seems good for logit models (Piironen and Vehtari, 2017a)
  real tau0 = (m0/(L-m0)) * (sig/sqrt(1.0*R*T)); // good if half of variables are expected to be zero
  real ss2 = 9; // slab scale squared
  real nu = 10; // determines degrees of freedom; implies marginally that beta ~ student-t(nu, 0, sqrt(ss2))
}

parameters {
  matrix[R, L] beta_tilde;
  matrix<lower=0>[R, L] lambda; // local prior for horseshoe
  real<lower=0> tau_tilde; // global prior for horseshoe
  real<lower=0> c2_tilde;
}

transformed parameters {
  matrix[R, L] B; // the matrix of coefficients
  {
    real tau = tau0 * tau_tilde; // tau ~ cauchy(0, tau0)
    real c2 = ss2 * c2_tilde; // c2 ~ inv_gamma(.5*nu, .5*nu*ss2)
    matrix[R, L] numer = (c2 * square(lambda));
    matrix[R, L] denom = (c2 + square(tau) * square(lambda));
    matrix[R, L] lambda_tilde = sqrt(numer ./ denom); // part of the horseshoe scale

    B = tau * lambda_tilde .* beta_tilde; // B ~ normal(0, tau*lambda_tilde)
  }
}

model {
  //priors
  to_vector(beta_tilde) ~ normal(0, 1);
  to_vector(lambda) ~ cauchy(0, 1);
  tau_tilde ~ cauchy(0, 1);
  c2_tilde ~ inv_gamma(0.5*nu, 0.5*nu);

  // likelihood
  for (r in 1:R) {
    for (t in 1:T) {
      Y[r, t] ~ categorical_logit(X[r,t]*B[r]');
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
