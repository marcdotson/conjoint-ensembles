// HBMNL for ANA pathology
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
  matrix[R, L] beta_tilde;
  matrix<lower=0>[R, L] lambda; // local prior for horseshoe
  real<lower=0> tau_tilde; // global prior for horseshoe
  real<lower=0> B_sig; // scale for B
}

transformed parameters {
  matrix[R, L] B = beta_tilde .* lambda * B_sig * tau_tilde;
}

model {
  //priors
  to_vector(beta_tilde) ~ normal(0, 1);
  to_vector(lambda) ~ cauchy(0, 1);
  tau_tilde ~ cauchy(0, 1);
  B_sig ~ normal(0, 2);

  // model fitting
  for (r in 1:R) {
    for (t in 1:T) {
      Y[r, t] ~ categorical_logit(X[r,t]*B[r]');
    }
  }
}

generated quantities {
  // Yp is predicted choices for new data.
  real Y_ppc[R, T];
  for (r in 1:R) {
    for (t in 1:T) {
      Y_ppc[r,t] = categorical_logit_rng(X[r,t] * B[r]');
    }
  }
}
