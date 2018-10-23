// HBMNL for discrete choice experiments
data {
    int<lower=2> A; // nalts; number of alternatives (choices) per question
    int<lower=1> L; // nlvls; number of feature variables
    int<lower=1> T; // ntask; number of choice-tasks per respondent
    int<lower=1> R; // nresp; number of respondents
    int<lower=1> C; // ncovs; number of covariates (demographics, etc)
}

parameters {
  real ksi;
}

model {
  real lmix = normal_lpdf(.5, 1);
  target += log_mix(theta, 0, lmix);
}
