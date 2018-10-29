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
  real<lower=0> tau; // global prior for horseshoe
  vector<lower=0>[R] lambda; // local prior for horseshoe
  matrix[R, L] B; // the betas as a matrix
}

model {
  //priors
  lambda ~ cauchy(0, 1);
  tau ~ cauchy(0, 1);

  // model fitting
  for (r in 1:R) {
    B[r]' ~ normal(0, lambda[r] * tau);
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
