// Finnish Horseshoe prior for multinomial logit
// Adapted from Michael Betancourt's case study on Bayes Sparse Regression
data {
  int<lower=2> A; // number of alternatives (choices) per question
  int<lower=1> L; // number of feature variables
  int<lower=1> N; // number of respondents
  int<lower=1> Ntest; // number of test cases

  matrix[A, L] X[N]; // matrix of attributes for each obs
  int<lower=1,upper=A> Y[N]; // observed responses

  matrix[A, L] Xtest[Ntest]; // test design matrix
}

transformed data {
  real m0 = 0.5*L;
  real sig = 2; // seems good for logit models (Piironen and Vehtari, 2017a)
  real tau0 = (m0/(L-m0)) * (sig/sqrt(1.0*N));
  real ss2 = 9; // slab scale squared
  real nu = 10; // determines degrees of freedom; implies marginally that beta ~ student-t(nu, 0, sqrt(ss2))
}

parameters {
  vector[L] beta_tilde;
  vector<lower=0>[L] lambda; // local prior for horseshoe
  real<lower=0> tau_tilde; // global prior for horseshoe
  real<lower=0> c2_tilde;
}

transformed parameters {
  vector[L] B; // the matrix of coefficients
  {
    real tau = tau0 * tau_tilde; // tau ~ cauchy(0, tau0)
    real c2 = ss2 * c2_tilde; // c2 ~ inv_gamma(.5*nu, .5*nu*ss2)
    vector[L] numer = (c2 * square(lambda));
    vector[L] denom = (c2 + square(tau) * square(lambda));
    vector[L] lambda_tilde = sqrt(numer ./ denom); // part of the horseshoe scale

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
  for (n in 1:N) {
    Y[n] ~ categorical_logit(X[n]*B);
  }
}

generated quantities {
  // Yp is predicted choices for new data.
  int<lower=0, upper=A> Yhat[Ntest];
  int<lower=0, upper=1> Yc[Ntest, A];

  for (n in 1:Ntest) {
    Yhat[n] = categorical_logit_rng(Xtest[n]*B);
    for (a in 1:A) {
      if (Yhat[n] == a) {
        Yc[n, a] = 1;
      }
      else {
        Yc[n, a] = 0;
      }
    }
  }
}

