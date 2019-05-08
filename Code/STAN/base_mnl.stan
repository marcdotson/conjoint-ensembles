// multinomial logit for discrete choice experiments
data {
  int<lower=2> A; // number of alternatives (choices) per question
  int<lower=1> L; // number of feature variables
  int<lower=1> R; // number of respondents
  int<lower=1> T; // number of tasks
  int<lower=1> Rtest; // number of test respondents
  int<lower=1> Ttest; // number of test tasks

  matrix[A, L] X[R, T]; // matrix of attributes for each obs
  int<lower=1, upper=A> Y[R, T]; // observed responses

  matrix[A, L] Xtest[Rtest, Ttest];

  real loc; // location of B
  real<lower=0> scale; // variance of B

}

transformed data {
  int<lower=1> N = R*T; // number of observations
  int<lower=1> Ntest = Rtest * Ttest; // number of test observations
  matrix[A, L] Xx[R*T]; // new x matrix
  int<lower=1, upper=A> Yy[R*T]; // new y matrix
  matrix[A, L] Xxtest[R*T]; // new test x matrix

  for (t in 1:T) {
    for (r in 1:R) {
      Xx[T*(r-1) + t] = X[r, t];
      Yy[T*(r-1) + t] = Y[r, t];
    }
  }

  for (t in 1:Ttest) {
    for (r in 1:Rtest) {
      Xxtest[Ttest*(r-1) + t] = Xtest[r, t];
    }
  }
}

parameters {
  vector[L] B; // matrix of beta coefficients
}

model {
  B ~ normal(loc, scale); // prior for B
  for (n in 1:N) {
    Yy[n] ~ categorical_logit(Xx[n]*B);
  }
}

generated quantities {
  // Yp is predicted choices for new data.
  int<lower=0, upper=A> Yhat[Ntest];
  int<lower=0, upper=1> Yc[Ntest, A];

  for (n in 1:Ntest) {
    Yhat[n] = categorical_logit_rng(Xxtest[n]*B);
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
