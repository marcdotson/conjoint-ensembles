// multinomial logit for discrete choice experiments
data {
  int<lower=2> A; // number of alternatives (choices) per question
  int<lower=1> L; // number of feature variables
  int<lower=1> N; // number of observations
  int<lower=1> N_test; // number of test observations
  int<lower=1, upper=A> Y[N]; // observed responses
  matrix[A, L] X[N]; // matrix of attributes for each obs
  real loc; // location of B
  real<lower=0> scale; // variance of B
  matrix[A, L] Xtest[N_test];
}

parameters {
  vector[L] B; // matrix of beta coefficients
}

model {
  B ~ normal(loc, scale); // prior for B
  for (n in 1:N) {
    Y[n] ~ categorical_logit(X[n]*B);
  }
}

generated quantities {
  // Yp is predicted choices for new data.
  int<lower=0, upper=A> Y_pp[N_test];
  int<lower=0, upper=1> Y_count[N_test, A];

  for (n in 1:N_test) {
    Y_pp[n] = categorical_logit_rng(Xtest[n]*B);
    for (a in 1:A) {
      if (Y_pp[n] == a) {
        Y_count[n, a] = 1;
      }
      else {
        Y_count[n, a] = 0;
      }
    }
  }
}
