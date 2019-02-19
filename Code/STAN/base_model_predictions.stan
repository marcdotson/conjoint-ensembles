// Generates predictions for a base model
data {
  int<lower=2> A; // number of alternatives (choices) per question
  int<lower=1> L; // number of feature variables
  int<lower=1> I; // number of iterations from the model
  int<lower=1> R; // number of respondents in training set
  int<lower=1> Rtest; // number of repondents in test set
  int<lower=1> Ttest; // number of choice tasks in test set
  real B[I, R, L]; // matrix of model coefficients
  matrix[A, L] Xtest[Rtest, Ttest]; // matrix of attributes for each obs
}

transformed data {
  vector[L] B_avg;

  for (l in 1:L) {
    B_avg[l] = mean(to_matrix(B[,,l]));
  }
}

generated quantities {
  int<lower=1, upper=A> Y_pp[Rtest, Ttest]; // predicted response
  int<lower=0, upper=1> Y_count[Rtest, Ttest, A]; // counts of response

  for (r in 1:Rtest) {
    for (t in 1:Ttest) {
      Y_pp[r, t] = categorical_logit_rng(Xtest[r, t] * B_avg);
      for (a in 1:A) {
        if (Y_pp[r, t] == a) {
          Y_count[r, t, a] = 1;
        }
        else {
          Y_count[r, t, a] = 0;
        }
      }
    }
  }
}
