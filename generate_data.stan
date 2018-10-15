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

    // unobservables
    matrix[R, L] Mu;
    matrix[L, R] alpha;
    matrix[L, L] L_Omega = lkj_corr_cholesky_rng(L, 5);
    vector<lower=0>[L] tau;
    matrix[R, L] Beta;

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

    Beta = Mu + (diag_pre_multiply(tau, L_Omega) * alpha)';

    for (r in 1:R) {
        for (t in 1:T) {
            // respondent r and task t chooses between random alternatives
            {
                vector[A*L] X_task;
                for (al in 1:A*L) {
                    X_task[al] = bernoulli_rng(.5);
                }
                X[r, t] = to_matrix(X_task, A, L);
            }
            Y[r, t] = categorical_logit_rng(X[r, t]*Beta[r]');
        }
    }
}
