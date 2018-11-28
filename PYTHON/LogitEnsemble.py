import numpy as np
import utils
from constants import *

data_dict = utils.get_data_dict()

stan_data = {
        'A':nalts,
        'L':nlvls,
        'loc':0,
        'scale':2
        }

step = np.linspace(0,N,ntask_train+1).astype(np.int64)

model_list = ['mnl', 'mnl_fhorseshoe']

M = len(model_list)

# fit model to data
Yhat_train = np.zeros((N,nalts,M))
for k in range(K):
    for m in range(M):
    
        k_fold = np.array([True]*N)
        k_fold[step[k]:step[k+1]] = False

        stan_data['N'] = sum(k_fold) # N-Nk
        stan_data['Ntest'] = sum(~k_fold) # Nk

        # new training set = subset of full training set
        stan_data['X'] = data_dict['Xtrain'][k_fold, :, :]
        stan_data['Y'] = data_dict['Ytrain'][k_fold]
    
        # new test set = complement of new training set | full training set
        stan_data['Xtest'] = data_dict['Xtrain'][~k_fold, :, :]
    
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
    stan_data['X'] = data_dict['Xtrain']
    stan_data['Y'] = data_dict['Ytrain']
 
    # full test set
    stan_data['Xtest'] = data_dict['Xtest']

    base_model = utils.get_model(model_name=model_list[m])
    FIT = utils.fit_model_to_data(base_model, stan_data)

    Yc_test = FIT.extract(pars=['Yc'])['Yc'].sum(axis=0)
    Yhat_k = np.argmax(Yc_test, axis=1)
    Yhat_test[np.array([True]*Ntest), Yhat_k, m] += 1

    model_scores.append(Ntest - np.count_nonzero(Yhat_k+1 - data_dict['Ytest']))


# Fit stacking model to full test data using augmented training set
stan_data['M'] = M
stan_data['Yhat_train'] = Yhat_train
stan_data['Yhat_test'] = Yhat_test

meta_model = utils.get_model(model_name='meta_mnl')
FIT = utils.fit_model_to_data(meta_model, stan_data)

Yc_stacking = FIT.extract(pars=['Yc'])['Yc'].sum(axis=0)
Yhat_stacking = np.argmax(Yc_stacking, axis=1) + 1

hit_count = Ntest - np.count_nonzero(Yhat_stacking - data_dict['Ytest'])
print("\n\nMODEL SCORES\n\t", np.array(model_scores)/Ntest)
print("ENSEMBLE SCORE\n\t",hit_count/Ntest)
print(Yhat_test)

## END ##
