import numpy as np
import utils

def get_data(choice=None):
    """
    """
    option_dict = {"01":"01_PathologyNone",
                   "02":"02_PathologyANA",
                   "03":"03_PathologyScreening",
                   "04":"04_PathologyMultiple"}

    if choice:
        try:
            return np.load("../../Data/{0}/data_dict.npz".format(option_dict[choice]))
        except:
            raise ValueError("Could not find dataset for string key {0}".format(choice))
    else:
        print("\n\nAvailable Datasets:\n\n")
        for d in sorted(option_dict.values()):
            print("\t",d)
        choice = input("\n\nenter a number as shown above>> ")
        assert choice in option_dict.keys()
    
        return np.load("../../Data/{0}/data_dict.npz".format(option_dict[choice]))


def hbmnl(data_dict, mu=(0,1), alpha=(0,10), lkj_param=5):
    """
    Hierarchical Bayesian Multi-Nomial Logit for conjoint analysis.

    INPUT
        data_dict (dict)

    OUTPUT
        results (dict)

    """

    # define local variables
    nresp = data_dict['X'].shape[0]
    nalts = data_dict['X'].shape[2]
    nlvls = data_dict['X'].shape[3]
    ntask_train = data_dict['Xtrain'].shape[1]
    ntask_test = data_dict['Xtest'].shape[1]
    N = nresp*ntask_train
    Ntest = nresp*ntask_test

    stan_data = {
            'A':nalts,
            'L':nlvls,
            'T':ntask_train,
            'R':nresp,
            'C':1,
            'Rtest':nresp,
            'Ttest':ntask_test,
            'X':data_dict['Xtrain'],
            'Y':data_dict['Ytrain'].astype(np.int64),
            'Z':np.ones((nresp,1)),
            'Xtest': data_dict['Xtest'],
            'mu_mean': mu[0],
            'mu_scale': mu[1],
            'alpha_mean': alpha[0],
            'alpha_scale': alpha[1],
            'lkj_param': lkj_param
            }
    
    
    # fit model to data
    base_model = utils.get_model(model_name='mnl_vanilla')
    FIT = utils.fit_model_to_data(base_model, stan_data)
    
    Yc = FIT.extract(pars=['Yc'])['Yc'].sum(axis=0).reshape((Ntest, nalts))
    Yhat = np.argmax(Yc, axis=1) + 1
    
    hit_count = Ntest - np.count_nonzero(Yhat - data_dict['Ytest'].reshape(Ntest))

    # store results
    results = dict()
    results["SCORE"] = hit_count/Ntest

    return results




def ensemble(data_dict):
    """
    Stacking Ensemble for conjoint analysis.

    INPUT
        data_dict (dict)

    OUTPUT
        results (dict)

    """

    # define local variables
    nresp = data_dict['X'].shape[0]
    ntask = data_dict['X'].shape[1]
    nalts = data_dict['X'].shape[2]
    nlvls = data_dict['X'].shape[3]
    N = nresp*data_dict['Xtrain'].shape[1]
    Ntest = nresp*data_dict['Xtest'].shape[1]
    K = 2
    M = nlvls
    
    # initialize and format data for stan model
    data = {}
    data['Xtrain'] = data_dict['Xtrain'].reshape(N, nalts, nlvls)
    data['Ytrain'] = data_dict['Ytrain'].reshape(N)
    data['Xtest'] = data_dict['Xtest'].reshape(Ntest, nalts, nlvls)
    data['Ytest'] = data_dict['Ytest'].reshape(Ntest)
    data['Ytrain'] = data['Ytrain'].astype(np.int64)
    data['Ytest'] = data['Ytest'].astype(np.int64)
    stan_data = {
            'A':nalts,
            'L':nlvls-1,
            'loc':0,
            'scale':2
            }
    
    # fit model to data
    step = [0, N//2, N]
    Yhat_train = np.zeros((N,nalts,M))
    for k in range(K):
        for m in range(M):
    
            k_fold = np.array([True]*N)
            k_fold[step[k]:step[k+1]] = False
    
            # set the K-fold temporary values of N and Ntest
            stan_data['N'] = sum(k_fold)
            stan_data['Ntest'] = sum(~k_fold)
    
            # new training set = subset of full training set
            # LOVO = Leave One Variable Out
            Xtrain_lovo = np.delete(data['Xtrain'][k_fold, :, :], m, 2)
            Xtest_lovo = np.delete(data['Xtrain'][~k_fold, :, :], m, 2)
    
            stan_data['X'] = Xtrain_lovo
            stan_data['Y'] = data['Ytrain'][k_fold]
        
            # new test set = complement of new training set | full training set
            stan_data['Xtest'] = Xtest_lovo
        
            base_model = utils.get_model(model_name='mnl')
            FIT = utils.fit_model_to_data(base_model, stan_data)
        
            Yc = FIT.extract(pars=['Yc'])['Yc'].sum(axis=0)
            Yhat_k = np.argmax(Yc, axis=1)
            Yhat_train[~k_fold, Yhat_k, m] += 1
    
    # make predictions on full test set using full training set
    model_scores = []
    Yhat_test = np.zeros((Ntest,nalts,M))
    for m in range(M):
    
        stan_data['N'] = N # length of training set
        stan_data['Ntest'] = Ntest # length of test set
        
        # full training set
        Xtrain_lovo = np.delete(data['Xtrain'], m, 2)
        Xtest_lovo = np.delete(data['Xtest'], m, 2)
    
        stan_data['X'] = Xtrain_lovo
        stan_data['Y'] = data['Ytrain']
     
        # full test set
        stan_data['Xtest'] = Xtest_lovo
    
        base_model = utils.get_model(model_name='mnl')
        FIT = utils.fit_model_to_data(base_model, stan_data)
    
        Yc_test = FIT.extract(pars=['Yc'])['Yc'].sum(axis=0)
        Yhat_k = np.argmax(Yc_test, axis=1)
        Yhat_test[np.array([True]*Ntest), Yhat_k, m] += 1
    
        model_scores.append(Ntest - np.count_nonzero(Yhat_k+1 - data['Ytest']))
    
    
    # Fit stacking model to full test data using augmented training set
    stan_data['M'] = M
    stan_data['Yhat_train'] = Yhat_train.copy()
    stan_data['Yhat_test'] = Yhat_test.copy()
    stan_data['L'] = nlvls
    
    meta_model = utils.get_model(model_name='meta_mnl')
    FIT = utils.fit_model_to_data(meta_model, stan_data)
    
    Yc_stacking = FIT.extract(pars=['Yc'])['Yc'].sum(axis=0)
    Yhat_stacking = np.argmax(Yc_stacking, axis=1) + 1
    model_weights = FIT.extract(pars=['B'])['B']
    
    ensemble_hit_count = Ntest - np.count_nonzero(Yhat_stacking - data['Ytest'])
    
    # store results
    results = dict()
    results["SCORE"] = ensemble_hit_count/Ntest
    results["BASE MODEL SCORES"] = np.array(model_scores)/Ntest
    results["MODEL WEIGHTS"] = np.around(model_weights.mean(axis=0),decimals=2)
    
    yy = Yhat_test.sum(axis=2)
    
    coverage_list = []
    for j in range(Ntest):
        coverage_list.append(max(yy[j, :]))
    coverage = np.array(coverage_list)
    
    model_coverage = np.zeros((M,2))
    for i in range(M):
        model_coverage[i,0] = i+1
        model_coverage[i,1] = len(coverage[coverage==i+1])
    results["BASE MODEL COVERAGE"] = model_coverage

    return results


## END ##
