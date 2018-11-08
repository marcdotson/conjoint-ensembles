data {
  int<lower=0> M; // number of models in ensemble
  int<lower=2> A; // number of alternatives
  int<lower=0> R; // number of test respondents
  int<lower=0> T; // number of test choice tasks
  vector<lower=0, upper=1>[M] w; // stacking weights
  vector<lower=1, upper=A>[A] Y_pp[M, R, T]; // posterior predictive density for model M
}

generated quantities {
  int Y_hat[R, T];
  for (t in 1:T) {
    for (r in 1:R) {
      Y_hat[r, t] = categorical_logit_rng(w .* Y_pp[, r, t]);
    }
  }
}
