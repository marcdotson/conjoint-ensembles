Code
================

## R

  - `01_simulate-data.R`: Simulate data with and without one or more
    pathologies present.
  - `02_clean-data.R`: Clean and restructure empirical application data
    for use in the ensemble.
  - `03_conjoint-ensemble.R`: Induce randomization, run the ensemble,
    and produce a consensus.
  - `04_competing-models.R`: Run the standard model or models specific
    to each pathology.
  - `05_model-comparison.R`: Compare fit across models.

## Python

  - `conjoint.py`: Functions for setting up data and running an HMNL,
    creating an ensemble, and stacking.
  - `demo.py`: Demo that calls on `conjoint.py`.
  - `test_utils.py`: Functions for testing.
  - `utils.py`: Utility functions.
  - `visualize.py`: Visualization functions.

## Stan

  - `hmnl_vanilla.stan`: Vanilla hierarchical multinomial logit with
    testing included in the `generated quantities` block.
  - `hmnl.stan`: Hierarchical multinomial logit.
  - `meta_mnl.stan`: Aggregate multinomial logit meta-learning model.
  - `mnl.stan`: Aggregate multinomial logit.
