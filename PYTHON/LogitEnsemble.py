import numpy as np
import utils
import matplotlib.pyplot as plt
from constants import *

# X has dimension R,T,A,L
# Y has dimension R,T

data_dict_hmnl,data_dict = utils.get_data_dict(kind='ensemble',pathology_type=pathology_type)

stan_data = {
        'A':nalts,
        'L':11,
        'loc':0,
        'scale':2
        }

step = np.linspace(0,N,K+1).astype(np.int64)

model_list = sorted(['mnl']*nlvls)

M = len(model_list)

# fit model to data
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
        Xtrain_lovo = np.delete(data_dict['Xtrain'][k_fold, :, :], m, 2)
        Xtest_lovo = np.delete(data_dict['Xtrain'][~k_fold, :, :], m, 2)

        stan_data['X'] = Xtrain_lovo
        stan_data['Y'] = data_dict['Ytrain'][k_fold]
    
        # new test set = complement of new training set | full training set
        stan_data['Xtest'] = Xtest_lovo
    
        base_model = utils.get_model(model_name=model_list[m])
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
    Xtrain_lovo = np.delete(data_dict['Xtrain'], m, 2)
    Xtest_lovo = np.delete(data_dict['Xtest'], m, 2)

    stan_data['X'] = Xtrain_lovo
    stan_data['Y'] = data_dict['Ytrain']
 
    # full test set
    stan_data['Xtest'] = Xtest_lovo

    base_model = utils.get_model(model_name=model_list[m])
    FIT = utils.fit_model_to_data(base_model, stan_data)

    Yc_test = FIT.extract(pars=['Yc'])['Yc'].sum(axis=0)
    Yhat_k = np.argmax(Yc_test, axis=1)
    Yhat_test[np.array([True]*Ntest), Yhat_k, m] += 1

    model_scores.append(Ntest - np.count_nonzero(Yhat_k+1 - data_dict['Ytest']))


# Fit stacking model to full test data using augmented training set
stan_data['M'] = M
stan_data['Yhat_train'] = Yhat_train.copy()
stan_data['Yhat_test'] = Yhat_test.copy()
stan_data['L'] = nlvls

meta_model = utils.get_model(model_name='meta_mnl2')
FIT = utils.fit_model_to_data(meta_model, stan_data)

Yc_stacking = FIT.extract(pars=['Yc'])['Yc'].sum(axis=0)
Yhat_stacking = np.argmax(Yc_stacking, axis=1) + 1
model_weights = FIT.extract(pars=['B'])['B']

ensemble_hit_count = Ntest - np.count_nonzero(Yhat_stacking - data_dict['Ytest'])

##################################
####### run standard model #######
##################################

stan_data = {
        'A':nalts,
        'L':nlvls,
        'T':ntask_train,
        'R':nresp_train,
        'C':ncovs,
        'Rtest':nresp_test,
        'Ttest':ntask_test,
        'X':data_dict_hmnl['Xtrain'],
        'Y':data_dict_hmnl['Ytrain'],
        'Z':data_dict_hmnl['Z'],
        'Xtest': data_dict_hmnl['Xtest']
        }


# fit model to data
base_model = utils.get_model(model_name='mnl_vanilla')
FIT = utils.fit_model_to_data(base_model, stan_data)


Yc = FIT.extract(pars=['Yc'])['Yc'].sum(axis=0).reshape((Ntest, nalts))
Yhat = np.argmax(Yc, axis=1) + 1

hit_count = Ntest - np.count_nonzero(Yhat - data_dict_hmnl['Ytest'].reshape(Ntest))
print("\nSTANDARD MODEL SCORE\n\t",hit_count/Ntest)
print("ENSEMBLE SCORE\n\t",ensemble_hit_count/Ntest)
print("BASE MODEL SCORES\n\t", np.array(model_scores)/Ntest)
print("MODEL WEIGHTS\n\t", np.around(model_weights.mean(axis=0),decimals=2))

yy = Yhat_test.sum(axis=2)

coverage_list = []
for j in range(Ntest):
    coverage_list.append(max(yy[j, :]))
coverage = np.array(coverage_list)

print("BASE MODEL COVERAGE")
for i in range(M):
    print(i+1, len(coverage[coverage==i+1]))

## END ##
