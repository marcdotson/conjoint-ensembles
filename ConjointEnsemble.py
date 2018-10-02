from PYTHON import generate_simulated_data as gsd
import matplotlib.pyplot as plt
import pystan
import pickle
import numpy as np


model_name = 'HBMNL_vanilla'

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
fit = sm.sampling(data=data_dict, iter=800, chains=2)
fit.extract(pars=['log_lik'])['log_lik']
