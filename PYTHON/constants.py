"""
Parameters and constant values for running HBMNL vs ensemble conjoint model.
"""

### SAMPLER PARAMETERS ###
niter = 300 # number of iterations in the HMC sampler (Stan param)
nchains = 2 # number of markov chains (mostly useful for diagnostics)
treedepth = 3 # how deep to explore the posterior space
random_seed = 1750532

### DATA SET PARAMETERS ###
nresp = 100 # number of respondents
ntask = 15 # number of tasks (combined training and holdout sets)
nalts = 4 # number of alternatives or options
nlvls = 12 # number of feature-levels
ncovs = 1 # number of covariates (for standard hbmnl)

nresp_train = 100 # number of respondents in the training set
ntask_train = 10 # number of choice-tasks in the training set
nresp_test = 100 # number of respondents in the test set
ntask_test = 5 # number of choice-tasks in the test set

pathology_type = 'basic' # describes the pathologies to induce in the data set
# Options:
#   'none' : don't induce pathologies on the data set
#   'screening' : induce a screening pathology (feature variable set to -inf)
#   'ANA' : induce attribute non-attendence (feature variable set to 0
#   'basic' : induce 'screening' and 'ANA'

K = 2 # K for k-fold CV

N = nresp_train*ntask_train # number of training observations
Ntest = nresp_train*ntask_test # number of test observations
