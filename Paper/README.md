Clever Randomization and Ensembling Strategies for Accommodating
Multiple Data Pathologies in Conjoint Studies
================

## Abstract

Respondent behavior in conjoint studies often deviates from the
assumptions of random utility theory. We refer to deviations from
normative choice behavior as data pathologies. A variety of models have
been developed that attempt to correct for specific pathologies (i.e.,
screening rules, respondent quality, attribute non-attendance, etc.).
While useful, these approaches tend to be both conceptually complex and
computationally intensive. As such, these approaches have not widely
diffused into the practice of marketing research. In this paper we draw
on innovations in machine learning to develop a practical approach that
relies on (clever) randomization strategies and ensembling to
simultaneously accommodate multiple data pathologies in a single model.
We provide tips and tricks on how to implement this approach in
practice.

Keywords: Choice Modeling, Machine Learning, Random Forests

License: This work is licensed under a [Creative Commons
Attribution-ShareAlike 4.0 International
License](https://creativecommons.org/licenses/by-sa/4.0/).

## Introduction

  - Define data pathologies
  - Prediction WRT conjoint

Ensemble-based approaches currently dominate the world of competitive
out-of-sample prediction. From Kaggle to the Netflix Prize, the
predictive power inherent in using many models overshadows prediction
reliant on the performance of a single model. The primary reason
ensembles predict so well is that they serve as a hedge against model
misspecification. Since we have uncertainty about the correct model for
any given context, running many models and producing a consensus is a
simple yet powerful way to improve predictions.

In the world of conjoint, most studies are conducted using a single
model. When the aim of a conjoint study is solely inference and not
prediction, a single-model approach is arguably best. The academic
literature for conjoint is filled with models designed to improve
inference, especially when respondents behave in ways that are
“pathological” to the standard model. However, there are three reasons
to argue for an ensemble-based approach to conjoint analysis. First, the
end goal of many conjoint studies is prediction in the form of accurate
market simulations. Second, we still have uncertainty about the correct
model for any given conjoint study. Third, there is no single model that
accounts for all the respondent behaviors that result in the “data
pathologies” that have been addressed separately in the literature.

The remainder of the paper will be organized as follows. In Section 2,
we walk through ensemble approaches to prediction. In Section 3, we
detail our ensemble approach to conjoint analysis. In Section 4, we
provide results from simulation studies and an empirical application. In
Section 5, we conclude.

### Previous Material

Many academic innovations in conjoint modeling in the recent past have
focused on developing models that capture individual data pathologies.
By data pathologies, we refer to respondent-specific behaviors that
deviate from one or more of the core assumptions of random utility
theory. Examples of these types of models include those that deal with
screening rules, attribute non-attendance, deviations from IIA choice
behavior, respondent quality, etc.

While each of these models have been shown to improve both model fit and
inference in contexts where a particular data pathology is present, none
of these approaches has had a meaningful impact on the practice of
conjoint. In our opinion, this is the result of three key factors.

  - The derivation and implementation of these models is complex. They
    are tough for practitioners to understand, implement and sell to
    clients.
  - There is a lack of commercial software than can be used to easily
    implement and simulate from these models. Custom coding is required
    that can dramatically increase the expense and cycle time of a
    project.  
  - Each model is designed to capture a single pathology in the data.
    While this is useful, it is highly likely that multiple pathologies
    may be present and problematic in any given study.

In this paper we explore alternative strategies for simultaneously
accommodating multiple data pathologies in a single modeling framework.
Our approach draws heavily on innovations in machine learning, in
particular the work of Breiman (2001) who showed remarkable improvement
in prediction resulting from averaging (or taking the consensus
prediction) from an ensemble of diverse models. Breiman’s approach
relies on the use of randomization to create diversity in the data used
to calibrate an ensemble of models, thus helping hedge against model
misspecification and over-fitting.

Our approach draws on the intuition of Brieman (2001) and is closely
related to Kevin Lattery’s 2015 Sawtooth Software Conference
presentation and corresponding paper. Lattery (2015) shows improvement
in predictive fit from an ensemble of diverse models that was created
using Latent Class methods. Our paper extends this work by considering
alternative strategies to create diversity in the ensemble of models
used for prediction. Specifically, we show that certain types of
randomization can be introduced into the construction of an ensemble in
order to mitigate particular data pathologies. For example, the type of
randomization used to accommodate the use of choice heuristics differs
from the type of randomization used to deal with respondent quality.

We will presents the theoretical justification for our approach and
practical strategies for implementation.

## Data Pathologies/Accommodating Data Pathologies

### Attribute Non-Attendance

### Screening Rules

## Ensemble Approaches to Prediction

Often, scientific models are used for explanation, meaning they are used
to test causal theories (Schmueli 2010). While important for social
science research, explanatory models may not yield good predictions.
Profit-maximizing firms depend on the ability to anticipate consumer
response to product offerings. Hence, what the firm truly desires is a
good predictive model. For conjoint analysis, the hierarchical linear
model is a good explanatory model but its predictions do not take full
advantage of the complexity in the data. This is particularly pronounced
when the data contains multiple data pathologies. Predictive accuracy
can be improved without totally sacrificing the ability to understand
why respondents made their choices.

One way to achieve greater predictive accuracy is to use a collection of
individual or base models aggregated in such a way that the aggregated
model performs better than any individual model with respect to some
performance criterion. This approach is called *ensembling* and was
formally introduced in (Breiman 2001). To see why an ensemble approach
is the preferable strategy for dealing with multiple data pathologies in
conjoint studies, consider the following alternative strategies.

The first alternative solution is to include the effects of the data
pathology in a specific model. Several papers have been published to
this effect \[*citation needed*\]. The advantage of this approach is
that the model will be the optimal model. However, this approach has
several disadvantages. The most obvious is that in order to account for
the pathology, the tailored model must increase in complexity. This
complexity is both technical and conceptual. More complex models are
less likely to be used by practitioners. Furthermore, this assumes that
the modeler is able to exactly specify the data generating process of
the pathology. This is a strong assumption and likely does not hold in
the majority of cases. In addition to the increased complexity, tailored
models are more computationally intensive than standard models. In this
case, the marginal benefit provided by the tailored model would have to
exceed the marginal cost associated with the additional computation.
Computational complexity usually increases super-linearly whereas
improved model performance is sub-linear or linear at best. Even if some
practitioners are able to account for the pathology and can afford the
increased computational cost, the tailored model is not robust to
pathologies in the data other than those accounted for by the model. If
other pathologies are unexpectedly present, the model will perform
suboptimally.

The second alternative solution is to account for data pathologies in a
way that is more flexible than the first solution. This could be
reasonably done during the process of constructing the priors of a
Bayesian model. This is the more principled modeling approach since
knowledge of potential data pathologies constitutes prior information
that should be represented by the prior distribution. The model can
represent any number of data pathologies as well as their interactions
by using a mixture prior. Additionally, we can still use the standard
multinomial logit model or any model of discrete choice that can be
formulated as a Bayesian model; accounting for pathologies in the prior
means we don’t have to change the structure of our model. With this, we
get accurate inference on all of the parameters jointly and the
resulting predictions are more reliable. Of course, the ideal case is
unrealistic and presents a variety of problems. Certainly, if only a few
pathologies are present, then the mixture model approach is best.
However, the computational intensity spikes for nontrivial pathologies.
In some cases, a sufficiently complicated prior may lead to an
intractable posterior, at least given the current tools for Bayesian
analysis. Furthermore, knowledge of the pathologies is expected to be
known beforehand. Sometimes this is reasonable but for general purpose
software solutions, this is unrealistic. In addition, the different
pathologies may interfere with each other in the prior. This can create
multi-modal posteriors (which results in model and computational
complexity) or unidentifiable regions in the posterior (which prevents
the model from capturing the effects of the pathology). Finally, adding
additional structure for pathologies that in reality are not present is
undesirable but largely unavoidable. This leaves us with a method that
is complete in scope but computationally unrealistic and does not scale
linearly with the addition of more pathologies.

Instead of the above alternatives, we advocate an ensemble approach. If
the ideal model is a Bayesian mixture model, we can approximate the
mixture model with an ensemble of base models fit to different subsets
of the input data. The critical assumption here is that the “true”
mixture model can be approximated as a convex combination of base
models. The goal is to choose a list of base models such that the ideal
mixture model is in the span of our list of base models. The way to
create this list is by assigning each base model to “cover” a specific
subset of the input data. This list of base models is chosen such that
their union forms a topological cover of the input space. Note that it
is not necessary to specify the exact base models that form the
theoretically true convex combination of base models; it is sufficient
to provide a list of models from which the convex combination can be
generated. Additionally, it is not necessary that the base models be
completely disjoint. However, the coverage of the ensemble of base
models will increase if the base model overlap is minimized.

By way of example, suppose each base model is trained on
![l-1](https://latex.codecogs.com/png.latex?l-1 "l-1") of the features
where ![l](https://latex.codecogs.com/png.latex?l "l") is total number
of independent variables in the data set. By restricting information,
each individual base model will under-perform relative to the standard
hierarchical model. However, this also allows the base model to gain
additional information by focusing more on the variables that are
available. When the base model predictions are aggregated and weighted
according to base model performance, the additional information
available will cause the ensemble to outperform the standard
hierarchical model in terms of out of sample predictive performance.

We could choose the base models a different way. For example, each base
models could be specified with a different prior that accounts for a
specific pathology or respondent behavior. Then the base models will
generate predictions that account for the effects of each pathology
respectively. Other ways of specifying base models include dropping
individuals or choice-tasks from the design matrix. In this paper, we
show results for the ensemble generated by base models trained on
![l-1](https://latex.codecogs.com/png.latex?l-1 "l-1") feature
variables.

The ensemble approach does come with some disadvantages. For example,
our inferences are certainly biased. However, a relatively small
increase in bias results in a much larger increase in predictive
accuracy (Schmueli 2010). The ensemble approach is also scalable in the
number of pathologies, easier to understand conceptually, and fairly
easy to implement. The ensemble will continue to improve as additional
base models are included with the existing base models.

## Model Specification

  - Describe Stacking as an ensemble approach
      - Point estimates versus distributions
  - Coverage metric as an approximate measure of base model overlap.
  - Describe how data was simulated and how models were compared

This paper applies the stacking approach to constructing an ensemble
predictor. Stacking is a method for combining the predictions of
![K](https://latex.codecogs.com/png.latex?K "K") different models by
weighting the models based on their performance in cross-validation
(Clarke 2018). In a stacking approach, a selection of base models is
chosen which are trained to input data
![X](https://latex.codecogs.com/png.latex?X "X") on outcomes
![Y](https://latex.codecogs.com/png.latex?Y "Y"). Each of these models
are trained on the input data and generate predictions for a holdout
set. The base models are tested for accuracy by comparing in-sample
predictions to the actual in-sample results. The base models are then
weighted according to their accuracy. Let
![X](https://latex.codecogs.com/png.latex?X "X") denote the experimental
design of the conjoint study and let
![Y](https://latex.codecogs.com/png.latex?Y "Y") denote the outcome. If
![f\_k](https://latex.codecogs.com/png.latex?f_k "f_k") is the
![k](https://latex.codecogs.com/png.latex?k "k")th base model with
corresponding parameters ![\\theta\_k \\in
\\Theta](https://latex.codecogs.com/png.latex?%5Ctheta_k%20%5Cin%20%5CTheta
"\\theta_k \\in \\Theta"), then the ensemble model
![F](https://latex.codecogs.com/png.latex?F "F") is

  
![F(X,Y|\\Theta) = \\sum\_{k=1}^K
\\hat{w\_k}f\_k(X,Y|\\theta\_k).](https://latex.codecogs.com/png.latex?F%28X%2CY%7C%5CTheta%29%20%3D%20%5Csum_%7Bk%3D1%7D%5EK%20%5Chat%7Bw_k%7Df_k%28X%2CY%7C%5Ctheta_k%29.
"F(X,Y|\\Theta) = \\sum_{k=1}^K \\hat{w_k}f_k(X,Y|\\theta_k).")  

The weights
![\\hat{w\_k}](https://latex.codecogs.com/png.latex?%5Chat%7Bw_k%7D
"\\hat{w_k}") can be optimized for predictive performance. according to
the logarithmic score criterion for the stacking of posterior
distributions given in (Yao et al 2018). Each of the
![f\_k](https://latex.codecogs.com/png.latex?f_k "f_k") are logit models
previously fit to a subset of the input data
![X,y](https://latex.codecogs.com/png.latex?X%2Cy "X,y"). We chose to
implement our model in Stan to take advantage of recent improvements to
the sampling strategy associated with Hamiltonian Monte Carlo. Details
on the model code can be found in the appendix.

  - Formal mathematical specification
      - Coverage Metric: The predictions generated by the ensemble are
        distributed among the possible alternatives for each choice
        task. Each of the models predicts which alternative is the one
        chosen by the respondent in the choice task. The coverage metric
        looks at the model predictions and for each observation
        (respondent-choice-task) and selects the maximum number of
        models whose predictions agree. This gives a rough sense of the
        heterogeneity in the model predictions. We want the models to be
        sufficiently similar that many agree but sufficiently different
        so that when they disagree, there is still a consensus. This is
        a rough measure of the topological cover argument we make above.

### Ensemble Procedure

Split data into training set, test set, and holdout set. The training
data will be used to train the base models. The holdout set is used to
determine ensemble performance (as opposed to the base model
performance). Divide the training set in half. Use the first half as a
simulated training set for the base models. Then generate predictions on
the input data of the second half of the training set. Repeat this but
switch the roles of the two halves of the training set.

Next, train the base models on the entire training set and make
predictions on the holdout set. You will now have base model predictions
for three data sets: the first half of the training set, the second half
of the training set, and the holdout set.

Concatenate the predictions from the two halves of the training set.
Regress (logistic) the actual responses on these base model predictions.
This gives coefficients that weight the base models according to their
ability to predict within the training set. Weight the base model
prediction distributions of the holdout set by the generated
coefficients. This generates new weighted prediction distributions. Find
the mode of this new predictive distribution. This mode is the
ensemble’s
prediction.

<img src="../Figures/Figure_1.png" title="The x-axis represents the induced pathology." alt="The x-axis represents the induced pathology." width="75%" style="display: block; margin: auto;" />

## Simulation Studies

<!-- \begin{table}[] -->

<!-- \centering -->

<!-- \caption{Model Performance Comparison} -->

<!-- \label{my-label} -->

<!-- \begin{tabular}{|l|l|l|l|l|l|} -->

<!-- \hline -->

<!-- \multicolumn{1}{|c|}{Cross-Validated Holdout Accuracy} & \multicolumn{4}{c|}{Data Generating Process} \\ \hline -->

<!-- \multicolumn{1}{|c|}{Model} & \multicolumn{1}{c|}{(0) Baseline} & \multicolumn{1}{c|}{(1) Screening} & \multicolumn{1}{c|}{(2) ANA} & \multicolumn{1}{c|}{(3) Interactions} \\ \hline -->

<!-- HB MNL & .7496 & .9252 & .5426 & .8342 \\ \hline -->

<!-- Screening &  &  &  &  \\ \hline -->

<!-- ANA &  &  &  &  \\ \hline -->

<!-- Proposed Ensemble & .7462 & .9190 & .5440 & .8398 \\ \hline -->

<!-- \end{tabular} -->

<!-- \caption{Based on the average of 10 runs, each run having a different data set. For a given run, the same data set was used for all models.} -->

<!-- \end{table} -->

## Empirical Application

## Conclusion

<!-- \newpage -->

<!-- \bibliographystyle{apalike} -->

<!-- %\bibliography{WIMI_Expedia} -->

<!-- %\appendix -->

<!-- %\newpage -->

<!-- %\section{Appendix: Brand Attitude Survey Questions}\label{app:survey_qs} -->

<!-- %\singlespacing -->

## Appendix

### Stan Model Code

``` stan
// HBMNL for discrete choice experiments
data {
  int<lower=2> A; // number of alternatives (choices) per question
  int<lower=1> L; // number of feature variables
  int<lower=1> R; // number of respondents
  int<lower=1> T; // number of questions (unique inquiries)
  int<lower=1> C; // number of respondent covariates (demographics, etc)
  int<lower=1> Rtest;
  int<lower=1> Ttest;

  matrix[A, L] X[R, T]; // matrix of attributes for each obs
  int<lower=1, upper=A> Y[R, T]; // observed responses
  matrix[R, C] Z; // vector of covariates for each respondent

  matrix[A, L] Xtest[Rtest, Ttest]; // test design matrix
}

parameters {
  matrix[L, R] alpha; // prior on variance of utilities B
  vector<lower=0, upper=pi()/2>[L] tau_unif;
  matrix[C, L] mu; // prior on mean of utilities B
  cholesky_factor_corr[L] L_Omega;

}

transformed parameters {
  matrix[R, L] B; // matrix of beta coefficients
  vector<lower=0>[L] tau; // prior scale
  for (l in 1:L) tau[l] = 2.5 * tan(tau_unif[l]);
  B = Z * mu + (diag_pre_multiply(tau,L_Omega) * alpha)';
}

model {
  //priors
  to_vector(alpha) ~ normal(0, 10);
  to_vector(mu) ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(5);


  // model fitting
  for (r in 1:R) {
    for (t in 1:T) {
      Y[r, t] ~ categorical_logit(X[r, t]*B[r]');
    }
  }
}

generated quantities {
  // Yp is predicted choices for new data.
  int<lower=0, upper=A> Yhat[Rtest, Ttest];
  int<lower=0, upper=1> Yc[Rtest, Ttest, A];
  vector[L] B_avg;

  for (r in 1:Rtest) {
    for (t in 1:Ttest) {
      for (l in 1:L) B_avg[l] = mean(B[,l]);
      Yhat[r, t] = categorical_logit_rng(Xtest[r, t]*B_avg);
      for (a in 1:A) {
        if (Yhat[r, t] == a) {
          Yc[r, t, a] = 1;
        }
        else {
          Yc[r, t, a] = 0;
        }
      }
    }
  }
}
```
