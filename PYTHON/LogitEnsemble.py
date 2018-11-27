import numpy as np
import utils
from constants import *

FIT = dict() # stores the separate model fits

data_dict = utils.generate_simulated_data(pathology_type='all')

data_dict['Xtrain'] = data_dict['X'][:nresp_train, :ntask_train, :, :].reshape(N, nalts, nlvls)
data_dict['Xtest'] = data_dict['X'][nresp_train:, ntask_train:, :, :].reshape(N_test, nalts, nlvls)

data_dict['Ytrain'] = data_dict['Y'][:nresp_train, :ntask_train].reshape(N)
data_dict['Ytest'] = data_dict['Y'][nresp_train:, ntask_train:].reshape(N_test)

stan_data = {
        'A':nalts,
        'L':nlvls,
        'N':N,
        'N_test':N_test,
        'X':data_dict['Xtrain'],
        'Y':data_dict['Ytrain'],
        'Xtest':data_dict['Xtest'],
        'Ytest':data_dict['Ytest'],
        'loc':0,
        'scale':2
        }


# fit model to data
Y_count_list = []
K = 10
for k in range(K):

    base_model = utils.get_model(model_name='mnl')
    FIT['{0}'.format(k)] = utils.fit_model_to_data(base_model, stan_data)

    Y_count_list.append(FIT['{0}'.format(k)].extract(pars=['Y_count'])['Y_count'].sum(axis=0))

print('')
for k in range(K):
    Y_count = Y_count_list[k]

    hit_count = 0
    for n in range(N_test):
        Y_predict = np.argmax(Y_count[n, :]) + 1
        if Y_predict == stan_data['Ytest'][n]:
            hit_count += 1
    
    print(k)
    print('\t', hit_count, '\t', hit_count/N_test, '\n')

Y_count = sum(Y_count_list)
hit_count = 0
for n in range(N_test):
    Y_predict = np.argmax(Y_count[n, :]) + 1
    if Y_predict == stan_data['Ytest'][n]:
        hit_count += 1

print('\n\nTOTAL')
print(hit_count, '/', N_test)
print(hit_count/N_test)

