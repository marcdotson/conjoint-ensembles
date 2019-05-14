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

parameters {
  vector[L] B; // matrix of beta coefficients
}

model {
  B ~ normal(loc, scale); // prior for B
  for (t in 1:T) {
    for (r in 1:R) {
      Y[r,t] ~ categorical_logit(X[r,t]*B);
    }
  }
}

generated quantities {
  // Yp is predicted choices for new data.
  int<lower=0, upper=A> Yhat[Rtest, Ttest];
  int<lower=0, upper=1> Yc[Rtest, Ttest, A];

  for (t in 1:Ttest) {
    for (r in 1:Rtest) {
      Yhat[r,t] = categorical_logit_rng(Xtest[r,t]*B);
      for (a in 1:A) {
        if (Yhat[r,t] == a) {
          Yc[r, t, a] = 1;
        }
        else {
          Yc[r, t, a] = 0;
        }
      }
    }
  }
}
