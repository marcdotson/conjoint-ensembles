Code
================

## R Scripts

The project work is broken down into the following steps with
corresponding scripts.

  - `01_simulate-data.R` Simulate data with and without pathologies and
    induce randomization.
  - `02_clean-data.R` Clean empirical application data and induce
    randomization.
  - `03_conjoint-ensemble.R` Run the ensemble using the clever
    randomization.
  - `04_meta-learner.R` Produce a consensus using the ensemble output.
  - `05_competing-models.R` Run the standard model or models specific to
    each pathology.
  - `06_model-comparison.R` Compare fit across models and ensembles.

## Functions

The following functions are stored in `Source`.

  - `ana_hmnl.R` Estimate an HMNL with attribute non-attendance.
  - `clever_randomization.R` Induce clever randomization and split into
    train and test data.
  - `conj_hmnl.R` Estimate an HMNL with conjunctive screening rules.
  - `conj_hmnl.cpp` Helper functions in C++ for `conj_hmnl.R`.
  - `ensemble_weights.R` Generate ensemble weights with a `stanfit`
    object as input.
  - `generate_data.stan` Generate data according to a HMNL.
  - `hmnl.R` Estimate an HMNL.
  - `hmnl.stan` Estimate an HMNL.
  - `hmnl_ensemble.stan` Estimate an HMNL as part of an ensemble.
  - `simulate_data.R` Simulate choice data with and without pathologies.
