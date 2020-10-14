Code
================

  - `01_simulate-data.R`: Simulate data with and without pathologies and
    induce randomization.
  - `02_clean-data.R`: Clean empirical application data and induce
    randomization.
  - `03_conjoint-ensemble.R`: Run the ensemble.
  - `04_meta-learner.R`: Produce a consensus using the ensemble output.
  - `05_competing-models.R`: Run the standard model or models specific
    to each pathology.
  - `06_model-comparison.R`: Compare fit across models and ensemble.

In `Source`:

  - `generate_data.stan`: Generate data according to the hierarchical
    multinomial logit.
  - `hmnl_noncentered.stan`: Hierarchical multinomial logit with a
    non-centered parameterization.
