import numpy as np
import utils
from constants import *

data_dict = utils.get_data_dict(kind='standard')

stan_data = {
        'A':nalts,
        'L':nlvls,
        'T':ntask_train,
        'R':nresp_train,
        'C':ncovs,
        'Rtest':nresp_test,
        'Ttest':ntask_test,
        'X':data_dict['Xtrain'],
        'Y':data_dict['Ytrain'],
        'Z':data_dict['Z'],
        'Xtest': data_dict['Xtest']
        }


# fit model to data
base_model = utils.get_model(model_name='mnl_vanilla')
FIT = utils.fit_model_to_data(base_model, stan_data)


Yc = FIT.extract(pars=['Yc'])['Yc'].sum(axis=0).reshape((Ntest, nalts))
Yhat = np.argmax(Yc, axis=1) + 1
print(Yhat)

hit_count = Ntest - np.count_nonzero(Yhat - data_dict['Ytest'].reshape(Ntest))
print("\nMODEL SCORE\n\t",hit_count/Ntest)
