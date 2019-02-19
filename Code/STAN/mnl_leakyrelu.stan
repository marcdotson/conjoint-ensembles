// HBMNL for discrete choice experiments
functions {
  vector LeakyReLU(vector X, int L) {
    vector[L] X_new;
    for (j in 1:L) {
      if (X[j] <= -1) X_new[j] = -square(X[j]);
      else X_new[j] = X[j];
    }
    return X_new;
  }
}

data {
  int<lower=2> A; // number of alternatives (choices) per question
  int<lower=1> L; // number of feature variables
  int<lower=1> N; // number of observations
  int<lower=1> Ntest; // number of test observations

  matrix[A, L] X[N]; // matrix of attributes for each obs
  int<lower=1,upper=A> Y[N]; // observed responses

  matrix[A, L] Xtest[Ntest]; // test design matrix
}

parameters {
  vector[L] beta; // coefficients
}

transformed parameters {
  vector[L] B; // matrix of beta coefficients
  B = LeakyReLU(beta, L);
}

model {
  //priors
  to_vector(beta) ~ normal(0, 1);

  // model fitting
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
