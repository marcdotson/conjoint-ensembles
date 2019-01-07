import numpy as np
import utils
from constants import *

# X is R,T,A,L
# Y is R,T

data_dict_hmnl,data_dict = utils.get_data_dict(kind='ensemble',pathology_type='basic')

stan_data = {
        'A':nalts,
        'L':nlvls//2,
        'loc':0,
        'scale':2
        }

step = np.linspace(0,N,ntask_train+1).astype(np.int64)
l_step = np.arange(0,nlvls+1,6)

model_list = sorted(['mnl', 'mnl_fhorseshoe']*2)

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
        stan_data['X'] = data_dict['Xtrain'][k_fold, :, l_step[m%2]:l_step[m%2 + 1]]
        stan_data['Y'] = data_dict['Ytrain'][k_fold]
    
        # new test set = complement of new training set | full training set
        stan_data['Xtest'] = data_dict['Xtrain'][~k_fold, :, l_step[m%2]:l_step[m%2 + 1]]
    
        base_model = utils.get_model(model_name=model_list[m])
        FIT = utils.fit_model_to_data(base_model, stan_data)
    
        Yc = FIT.extract(pars=['Yc'])['Yc'].sum(axis=0)
        Yhat_k = np.argmax(Yc, axis=1)
        Yhat_train[~k_fold, Yhat_k, m] += 1
        #Yhat_train[~k_fold, :, m] = Yc


# make predictions on full test set using full training set
model_scores = []
Yhat_test = np.zeros((Ntest,nalts,M))
for m in range(M):

    stan_data['N'] = N # length of training set
    stan_data['Ntest'] = Ntest # length of test set
    
    # full training set
    stan_data['X'] = data_dict['Xtrain'][:, :, l_step[m%2]:l_step[m%2 + 1]]
    stan_data['Y'] = data_dict['Ytrain']
 
    # full test set
    stan_data['Xtest'] = data_dict['Xtest'][:,:,l_step[m%2]:l_step[m%2 + 1]]

    base_model = utils.get_model(model_name=model_list[m])
    FIT = utils.fit_model_to_data(base_model, stan_data)

    Yc_test = FIT.extract(pars=['Yc'])['Yc'].sum(axis=0)
    Yhat_k = np.argmax(Yc_test, axis=1)
    Yhat_test[np.array([True]*Ntest), Yhat_k, m] += 1
    #Yhat_test[:, :, m] += Yc_test

    model_scores.append(Ntest - np.count_nonzero(Yhat_k+1 - data_dict['Ytest']))


# Fit stacking model to full test data using augmented training set
stan_data['M'] = M
stan_data['Yhat_train'] = Yhat_train
stan_data['Yhat_test'] = Yhat_test
stan_data['L'] = nlvls

meta_model = utils.get_model(model_name='meta_mnl2')
FIT = utils.fit_model_to_data(meta_model, stan_data)

Yc_stacking = FIT.extract(pars=['Yc'])['Yc'].sum(axis=0)
Yhat_stacking = np.argmax(Yc_stacking, axis=1) + 1

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
print(Yhat)

hit_count = Ntest - np.count_nonzero(Yhat - data_dict_hmnl['Ytest'].reshape(Ntest))
print("\nMODEL SCORE\n\t",hit_count/Ntest)
print("ENSEMBLE SCORE\n\t",ensemble_hit_count/Ntest)
print("BASE MODEL SCORES\n\t", np.array(model_scores)/Ntest)

yy = Yhat_test.sum(axis=2)

coverage_list = []
for j in range(Ntest):
    coverage_list.append(max(yy[j, :]))
coverage = np.array(coverage_list)

print("BASE MODEL COVERAGE")
for i in range(M):
    print(i+1, len(coverage[coverage==i+1]))
print(Yhat_test)

## END ##
