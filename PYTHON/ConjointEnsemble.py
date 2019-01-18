from sklearn.linear_model import LogisticRegression
import numpy as np
import utils
from constants import *

data_dict_hmnl,data_dict = utils.get_data_dict(kind='ensemble',pathology_type='none')

stan_data = {
        'A':nalts,
        'L':11,
        'loc':0,
        'scale':2
        }

step = np.linspace(0,N,ntask_train+1).astype(np.int64)

M = nlvls

# fit model to data
Yhat_train = np.zeros((N,nalts,M))
Yhats_train = np.zeros((N,M))
for k in range(K):
    for m in range(M):

        k_fold = np.array([True]*N)
        k_fold[step[k]:step[k+1]] = False

        stan_data['N'] = sum(k_fold) # N-Nk
        stan_data['Ntest'] = sum(~k_fold) # Nk

        # new training set = subset of full training set
        # LOVO = Leave One Variable Out
        Xtrain_lovo = np.delete(data_dict['Xtrain'][k_fold, :, :], m, 2)
        Xtest_lovo = np.delete(data_dict['Xtrain'][~k_fold, :, :], m, 2)

        stan_data['X'] = Xtrain_lovo.reshape((stan_data['N'],nalts*(nlvls-1)))
        stan_data['Y'] = data_dict['Ytrain'][k_fold].reshape(stan_data['N'])
    
        # new test set = complement of new training set | full training set
        stan_data['Xtest'] = Xtest_lovo.reshape((stan_data['Ntest'],nalts*(nlvls-1)))
        
        base_model = LogisticRegression(random_state=0, solver='lbfgs', multi_class='multinomial')
        FIT = base_model.fit(stan_data['X'], stan_data['Y'])
        Yhat_k = FIT.predict(stan_data['Xtest']) - 1

        #base_model = utils.get_model(model_name=model_list[m])
        #FIT = utils.fit_model_to_data(base_model, stan_data)
        #Yc = FIT.extract(pars=['Yc'])['Yc'].sum(axis=0)
        #Yhat_k = np.argmax(Yc, axis=1)

        Yhat_train[~k_fold, Yhat_k, m] += 1
        Yhats_train[~k_fold, m] += Yhat_k+1


# make predictions on full test set using full training set
model_scores = []
Yhat_test = np.zeros((Ntest,nalts,M))
Yhats_test = np.zeros((Ntest,M))
for m in range(M):

    stan_data['N'] = N # length of training set
    stan_data['Ntest'] = Ntest # length of test set
    
    # full training set
    Xtrain_lovo = np.delete(data_dict['Xtrain'], m, 2)
    Xtest_lovo = np.delete(data_dict['Xtest'], m, 2)

    stan_data['X'] = Xtrain_lovo.reshape((stan_data['N'],nalts*(nlvls-1)))
    stan_data['Y'] = data_dict['Ytrain'].reshape(stan_data['N'])
 
    # full test set
    stan_data['Xtest'] = Xtest_lovo.reshape((stan_data['Ntest'],nalts*(nlvls-1)))

    base_model = LogisticRegression(random_state=0, solver='lbfgs', multi_class='multinomial')
    FIT = base_model.fit(stan_data['X'], stan_data['Y'])
    Yhat_k = FIT.predict(stan_data['Xtest']) - 1


    #base_model = utils.get_model(model_name=model_list[m])
    #FIT = utils.fit_model_to_data(base_model, stan_data)
    #Yc_test = FIT.extract(pars=['Yc'])['Yc'].sum(axis=0)
    #Yhat_k = np.argmax(Yc_test, axis=1)

    Yhat_test[np.array([True]*Ntest), Yhat_k, m] += 1
    Yhats_test[:,m] += Yhat_k

    model_scores.append(Ntest - np.count_nonzero(Yhat_k+1 - data_dict['Ytest']))


# Fit stacking model to full test data using augmented training set

stan_data['M'] = M
stan_data['Yhat_train'] = Yhats_train
stan_data['Yhat_test'] = Yhats_test
stan_data['Y'] = data_dict['Ytrain'].flatten()
stan_data['L'] = nlvls


meta_model = LogisticRegression(random_state=0, solver='lbfgs', multi_class='multinomial')
FIT = meta_model.fit(stan_data['Yhat_train'], stan_data['Y'])
Yhat_stacking = FIT.predict(stan_data['Yhat_test'])

#print(Yhat_test.sum(axis=2))
#print(data_dict['Ytest'])
#input()

Phat = np.zeros((Yhat_test.shape[0],Yhat_test.shape[1]))
for m in range(M):
    Phat += Yhat_test[:,:,m]*model_scores[m]

Yhat_stacking = np.argmax(Phat, axis=1) + 1
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


base_model = utils.get_model(model_name='mnl_vanilla')
FIT = utils.fit_model_to_data(base_model, stan_data)


Yc = FIT.extract(pars=['Yc'])['Yc'].sum(axis=0).reshape((Ntest, nalts))
Yhat = np.argmax(Yc, axis=1) + 1

hit_count = Ntest - np.count_nonzero(Yhat - data_dict['Ytest'].reshape(Ntest))
print("\nSTANDARD MODEL SCORE\n\t",hit_count/Ntest)
print("ENSEMBLE SCORE\n\t", ensemble_hit_count/Ntest)
print("BASE MODEL SCORES\n\t", np.array(model_scores)/Ntest)

yy = Yhat_test.sum(axis=2)

coverage_list = []
for j in range(Ntest):
    coverage_list.append(max(yy[j, :]))
coverage = np.array(coverage_list)

print("BASE MODEL COVERAGE")
for i in range(M):
    print(i+1, len(coverage[coverage==i+1]))

## END ##
