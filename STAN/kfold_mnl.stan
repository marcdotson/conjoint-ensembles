// multinomial logit for discrete choice experiments
data {
  int<lower=2> A; // number of alternatives (choices) per question
  int<lower=1> L; // number of feature variables
  int<lower=1> N; // number of observations
  int<lower=1> Nk; // number of CV fold k observations
  int<lower=1> Ntest; // number of out of sample observations

  matrix[A, L] X[N]; // matrix of attributes for each obs
  int<lower=1, upper=A> Y[N]; // observed responses

  matrix[A, L] Xk[Nk];
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
  int<lower=0, upper=A> Yhat[Nk];
  int<lower=0, upper=1> Yc[Nk, A];
  int<lower=0, upper=A> Yhat_test[Ntest];
  int<lower=0, upper=1> Yc_test[Ntest, A];

  for (n in 1:Nk) {
    Yhat[n] = categorical_logit_rng(Xk[n]*B);
    for (a in 1:A) {
      if (Yhat[n] == a) {
        Yc[n, a] = 1;
      }
      else {
        Yc[n, a] = 0;
      }
    }
  }

  for (n in 1:Ntest) {
    Yhat_test[n] = categorical_logit_rng(Xtest[n]*B);
    for (a in 1:A) {
      if (Yhat_test[n] == a) {
        Yc_test[n, a] = 1;
      }
      else {
        Yc_test[n, a] = 0;
      }
    }
  }

}
