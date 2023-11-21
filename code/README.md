# Code

## R Scripts

The project work is broken down into the following steps with
corresponding scripts.

- `01_control-file.R` Control for running the subsequent individual
  scripts.
- `02_data-prep.R` Simulate data or clean empirical data and induce
  randomization.
- `03_conjoint-ensemble.R` Run the conjoint ensemble using the clever
  randomization.
- `04_meta-learner.R` Produce weights using the ensemble output.
- `05_competing-models.R` Run the models specific to the indicated
  pathology.
- `06_model-comparison.R` Compute and compare fit across models.

## Functions

The following functions are stored in `/code/src`.

- `ana_hbmnl.R` One of two functions to estimate an HMNL with attribute
  non-attendance.
- `ana_hmnl.R` One of two functions to estimate an HMNL with attribute
  non-attendance.
- `conj_hmnl.R` Estimate an HMNL with conjunctive screening rules.
- `conj_hmnl.cpp` Helper functions in C++ for `conj_hmnl.R`.
- `ensemble_weights.R` Generate ensemble weights with a `stanfit` object
  as input.
- `generate_data.stan` Generate data according to a HMNL.
- `hit_prob.R` Compute hit probabilities.
- `hit_rate.R` Compute hit rates.
- `hmnl_ensemble.stan` Estimate an HMNL as part of an ensemble.
- `hmnl.R` Estimate an HMNL.
- `hmnl.stan` Estimate an HMNL.
- `predictive_fit_ensemble.R` Compute predictive fit for an ensemble of
  models.
- `predictive_fit_hmnl.R` Compute predictive fit for an HMNL.
