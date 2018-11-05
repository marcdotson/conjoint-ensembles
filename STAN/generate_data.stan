// HBMNL for discrete choice experiments
data {
  int<lower=2> A; // nalts; number of alternatives (choices) per question
  int<lower=1> L; // nlvls; number of feature variables
  int<lower=1> T; // ntask; number of choice-tasks per respondent
  int<lower=1> R; // nresp; number of respondents
  int<lower=1> C; // ncovs; number of covariates (demographics, etc)
}

generated quantities {
  // observables
  int<lower=1, upper=A> Y[R, T]; // outcome of choice task T for respondent R
  matrix[A, L] X[R, T]; // design matrix
  matrix[R, C] Z; //covariates

  // unobservables
  matrix[R, L] Mu; // prior mean for respondent coefficients
  matrix[L, R] alpha; // prior variance vector for respondent coefficients
  matrix[L, L] L_Omega = lkj_corr_cholesky_rng(L, 2); // correlation matrix
  vector<lower=0>[L] tau; // parameter of the covariance matrix
  matrix[R, L] B; // respondent coefficients

  // temporarily fix demographics to be constant
  for (r in 1:R) Z[r] = rep_row_vector(1, C);

  {
    vector[R*L] Mu_temp;
    vector[L*R] alpha_temp;

    for (rl in 1:R*L) {
      Mu_temp[rl] = normal_rng(0,1);
      alpha_temp[rl] = normal_rng(0,10);
    }
    Mu = to_matrix(Mu_temp, R, L);
    alpha = to_matrix(alpha_temp, L, R);
  }

  for (l in 1:L) tau[l] = 2.5 * tan(beta_rng(1, 1));

  B = Mu + (diag_pre_multiply(tau, L_Omega) * alpha)';

  // iterate through each respondent
  for (r in 1:R) {
    // assign pathologies
    {
      // randomly assign attribute non-attendance pathology
      if (bernoulli_rng(.5) == 1) {
        for (l in 1:L) {
          int p = bernoulli_rng(.5);
          B[r, l] *= p;
        }
      }
      // randomly assign screening behavior pathology
      else if (bernoulli_rng(.5) == 1) {
        for (l in 1:L) {
          int p = bernoulli_rng(.5);
          B[r, l] -= 100*p;
        }
      }
    }

    // now that the respondent coefficients are set, iterate through the choice tasks
    for (t in 1:T) {
      {
        vector[A*L] X_task;

        // randomly assign feature-levels in design matrix
        for (al in 1:A*L) {
          X_task[al] = bernoulli_rng(.5);
        }

        // format design matrix to have dimensions AxL
        X[r, t] = to_matrix(X_task, A, L);
      }

      // on task t, respondent r chooses among the alternatives
      Y[r, t] = categorical_logit_rng(X[r, t]*B[r]');
    }
  }
}
