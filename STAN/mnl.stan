// multinomial logit for discrete choice experiments
data {
  int<lower=2> A; // number of alternatives (choices) per question
  int<lower=1> L; // number of feature variables
  int<lower=1> N; // number of observations
  int<lower=1> Ntest; // number of out of sample observations

  matrix[A, L] X[N]; // matrix of attributes for each obs
  int<lower=1, upper=A> Y[N]; // observed responses

  matrix[A, L] Xtest[Ntest];

  real loc; // location of B
  real<lower=0> scale; // variance of B

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
