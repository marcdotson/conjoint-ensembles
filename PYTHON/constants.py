nresp = 100
ntask = 15
nalts = 4
nlvls = 12
ncovs = 1
niter = 300
nchains = 2
treedepth = 3
random_seed = 1750532

nresp_train = 100
ntask_train = 10
#nresp_test = 100
ntask_test = 5

K = ntask_train # K for k-fold CV

N = nresp_train*ntask_train
Nk = nresp_train
Ntest = nresp_train*ntask_test
