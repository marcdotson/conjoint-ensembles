from PYTHON import utils
import pystan

def get_stan_model(fname):
    with open(fname, 'r') as f:
        stan_model = f.read()
    sm = pystan.StanModel(model_code=stan_model)
    return sm

# generate simulated data
#sm1 = get_stan_model("generate_data.stan")

#data_dict = {'A':4, 'L':12, 'T':10, 'R':100, 'C':1}
#fit1 = sm1.sampling(data=data_dict,iter=100,warmup=0,chains=1,refresh=100,algorithm="Fixed_param")
#params1 = [i for i in fit1.sim['samples'][0].chains.values()][:-1]
#params1 = fit1.extract(pars=['X','Y','Z'])
data_dict = utils.generate_simulated_data()

sm2 = get_stan_model("HBMNL_vanilla.stan")
fit2 = sm2.sampling(data=data_dict)
#params2 = [i for i in fit2.sim['samples'][0].chains.values()][:-1]

Y_ppc = fit2.extract(pars=['Y_ppc'])['Y_ppc']

