import numpy as np
import utils
from constants import *

FIT = dict() # stores the separate model fits

simulated_data = utils.generate_simulated_data(pathology_type='all')

train_data = {
        'A':nalts,
        'L':nlvls,
        'T':ntask_train,
        'R':nresp_train,
        'C':ncovs,
        'X':simulated_data['X'][:nresp_train, :ntask_train, :, :],
        'Y':simulated_data['Y'][:nresp_train, :ntask_train],
        'Z':simulated_data['Z'][:nresp_train, :]
        }


test_data = {
        'A':nalts,
        'L':nlvls,
        'I':int((nchains/2)*niter),
        'R':nresp_train,
        'Ttest':ntask_test,
        'Rtest':nresp_test,
        'C':ncovs,
        'Xtest':simulated_data['X'][nresp_train:, ntask_train:, :, :],
        'Y':simulated_data['Y'][nresp_train:, ntask_train:],
        'Z':simulated_data['Z'][nresp_train:, :]
        }


# fit model to data
base_model = utils.get_model(model_name='mnl_vanilla')
FIT['train'] = utils.fit_model_to_data(base_model, train_data)

test_data['B'] = FIT['train'].extract(pars=['B'])['B']
base_model_predict = utils.get_model(model_name='base_model_predictions')
FIT['test'] = utils.fit_model_to_data(base_model_predict, test_data, sampling_alg='Fixed_param')


Y_count_array = FIT['test'].extract(pars=['Y_count'])['Y_count']
Y_count = Y_count_array.sum(axis=0)

hit_count = 0

for r in range(nresp_test):
    for t in range(ntask_test):
        Y_predict = np.argmax(Y_count[r, t, :]) + 1
        if Y_predict == test_data['Y'][r, t]:
            hit_count += 1

print(hit_count, N_test)
print(hit_count/N_test)
