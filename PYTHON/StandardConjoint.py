import numpy as np
import utils
import matplotlib.pyplot as plt
from constants import *

# X has dimension R,T,A,L
# Y has dimension R,T

data_dict = utils.load_data_dict("./DATA/09_PathologyMultiple/data_dict.npz")


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
        'X':data_dict['Xtrain'],
        'Y':data_dict['Ytrain'],
        'Z':np.ones((nresp,1)),
        'Xtest': data_dict['Xtest']
        }


# fit model to data
base_model = utils.get_model(model_name='mnl_vanilla')
FIT = utils.fit_model_to_data(base_model, stan_data)


Yc = FIT.extract(pars=['Yc'])['Yc'].sum(axis=0).reshape((Ntest, nalts))
Yhat = np.argmax(Yc, axis=1) + 1

hit_count = Ntest - np.count_nonzero(Yhat - data_dict['Ytest'].reshape(Ntest))
print("\nSTANDARD MODEL SCORE\n\t",hit_count/Ntest)

## END ##
