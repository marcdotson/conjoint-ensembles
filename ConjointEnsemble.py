from PYTHON import generate_simulated_data as gsd
import matplotlib.pyplot as plt
import pystan
import pickle
import numpy as np


model_name = 'HBMNL_02'

with open('./MODELS/{0}.stan'.format(model_name), 'r') as f:
    stan_model = f.read()

try:
    sm = pickle.load(open('./MODELS/{0}.pkl'.format(model_name), 'rb'))

except:
    sm = pystan.StanModel(model_code=stan_model)
    with open('./MODELS/{0}.pkl'.format(model_name), 'wb') as f:
        pickle.dump(sm, f)

# ANA pathology
data_dict = gsd.generate_simulated_data(pathology_type='ANA')
data_dict['w'] = np.random.binomial(1, .5, size=data_dict['K'])
data_dict['Q'] = np.zeros((data_dict['J'], data_dict['K']))
fit = sm.sampling(data=data_dict, iter=800, chains=2)
fit.extract(pars=['log_lik'])['log_lik']


# Screening Pathology
#data_dict = gsd.generate_simulated_data(pathology_type='screening')
#data_dict['w'] = np.ones(data_dict['K'])
#data_dict['Q'] = np.random.choice([0,-np.inf], size=data_dict['J']*data_dict['K'], p=[.8, .2]).reshape((data_dict['J'], data_dict['K']))
#fit = sm.sampling(data=data_dict, iter=800, chains=2)
#fit.extract(pars=['log_lik'])['log_lik']
#print(fit)
