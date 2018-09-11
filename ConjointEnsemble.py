from PYTHON import generate_simulated_data as gsd
import matplotlib.pyplot as plt
import pystan
import pickle
import numpy as np

data_dict = gsd.generate_simulated_data()

model_name = 'HBMNL_01'

with open('./MODELS/{0}.stan'.format(model_name), 'r') as f:
    stan_model = f.read()

try:
    sm = pickle.load(open('./MODELS/{0}.pkl'.format(model_name), 'rb'))

except:
    sm = pystan.StanModel(model_code=stan_model)
    with open('./MODELS/{0}.pkl'.format(model_name), 'wb') as f:
        pickle.dump(sm, f)

fit = sm.sampling(data=data_dict, iter=800, chains=2)
fit.extract(pars=['log_lik'])['log_lik']
