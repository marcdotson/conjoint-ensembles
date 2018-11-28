// multinomial logit for discrete choice experiments
data {
  int<lower=1> M; // number of base models
  int<lower=2> A; // number of alternatives (choices) per question
  int<lower=1> L; // number of feature levels
  int<lower=1> N; // number of observations
  int<lower=1> Ntest; // number of out of sample observations

  matrix[A, L] X[N]; // matrix of attributes for each obs
  int<lower=1, upper=A> Y[N]; // observed responses

  matrix[A, L] Xtest[Ntest]; // test design matrix
  matrix[A, M] Yhat_train[N]; // base model predictions (or distributions)
  matrix[A, M] Yhat_test[Ntest]; // base model predictions (or distributions)

  real loc; // location of B
  real<lower=0> scale; // variance of B
}

transformed data {
  matrix[A, L+M] XYhat_train[N];
  matrix[A, L+M] XYhat_test[Ntest];
  for (n in 1:N) {
    XYhat_train[n] = append_col(X[n], Yhat_train[n]);
  }
  for (n in 1:Ntest) {
    XYhat_test[n] = append_col(Xtest[n], Yhat_test[n]);
  }
}

parameters {
  vector[L+M] B; // matrix of beta coefficients
}

model {
  B ~ normal(loc, scale); // prior for B
  for (n in 1:N) {
    Y[n] ~ categorical_logit(XYhat_train[n]*B);
  }
}

generated quantities {
  // Yhat is predicted choices for new data.
  int<lower=0, upper=A> Yhat[Ntest];
  int<lower=0, upper=1> Yc[Ntest, A];

  for (n in 1:Ntest) {
    Yhat[n] = categorical_logit_rng(XYhat_test[n]*B);
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
