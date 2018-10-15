import pystan

data_dict = {'A':4, 'L':12, 'T':10, 'R':100, 'C':1}
with open("generate_data.stan", 'r') as f:
    stan_model = f.read()

sm = pystan.StanModel(model_code=stan_model)
fit = sm.sampling(data=data_dict,iter=100,warmup=0,chains=1,refresh=100,algorithm="Fixed_param")
