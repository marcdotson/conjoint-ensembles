import pystan
import pickle
import numpy as np


### STAN SAMPLER PARAMETERS ###
niter = 300 # number of iterations in the HMC sampler (Stan param)
nchains = 2 # number of markov chains (mostly useful for diagnostics)
treedepth = 3 # how deep to explore the posterior space
random_seed = 1750532
np.random.seed(seed=random_seed)

    
def pathology(beta, kind="none", prob=[.5, .5]):
    # apply the pathologies
    if kind == 'none':
        pass
    elif kind == 'ANA':
        beta *= np.random.choice([1, 0], size=len(beta), p=prob)
    elif kind == 'screening':
        beta[-3:] -= 100
    elif kind == 'exponential':
        beta = np.random.exponential(size=len(beta))
    elif kind == 'screening_random':
        if int(np.random.choice([1,0], p=prob)):
            i = np.random.randint(len(beta))
            beta[i] = -10000
    elif kind == 'ANA_systematic':
        pathology_vector = np.ones_like(beta)
        pathology_vector[int(len(beta)//2):] = 0
        beta *= pathology_vector
    elif kind == 'ANA_random':
        if int(np.random.choice([1, 0], p=prob)):
            pathology_vector = np.ones_like(beta)
            pathology_vector[:int(len(beta)//2)] = 0
            beta *= pathology_vector
    elif kind == 'basic':
        beta = pathology(beta, kind='ANA')
        beta = pathology(beta, kind='screening')
    elif kind == 'all':
        # ANA
        if int(np.random.choice([1,0], p=prob)):
            beta *= np.random.choice([1, 0], size=len(beta), p=prob)
        # Screening
        if int(np.random.choice([0,1], p=prob)):
            beta[-3:] -= 100
        if int(np.random.choice([0,1], p=prob)):
            beta *= np.random.uniform(-10, 10, size=len(beta))
        if int(np.random.choice([0,1], p=prob)):
            beta = np.random.laplace(size=len(beta))
        if int(np.random.choice([0,1], p=prob)):
            beta += np.random.exponential(size=len(beta))

    return beta


def generate_simulated_design():
    # X is the experimental design
    X = np.zeros((nresp, ntask, nalts, nlvls))
    # Z is a matrix for demographic attributes
    Z = np.zeros((ncovs, nresp))
    
    for resp in range(nresp):
        z_resp = 1
        if ncovs > 1:
            raise NotImplementedError
    
        for scn in range(ntask):
            X_scn = np.random.choice([0,1], p =[.5,.5], size=nalts*nlvls).reshape(nalts,nlvls)
            X[resp, scn] += X_scn
    
        Z[:, resp] += z_resp
    
    # dictionary to store the simulated data and generation parameters
    data_dict = {'X':X,
                 'Z':Z.T,
                 'A':nalts,
                 'R':nresp,
                 'C':ncovs,
                 'T':ntask,
                 'L':nlvls}
    return data_dict


def compute_beta_response(data_dict, pathology_type=None):

    # beta means
    Gamma = np.random.uniform(-3, 4, size=data_dict['C'] * data_dict['L'])
    # beta variance-covariance
    Vbeta = np.diag(np.ones(data_dict['L'])) + .5 * np.ones((data_dict['L'], data_dict['L']))

    # Y is the response
    Y = np.zeros((data_dict['R'], data_dict['T']))
    # Beta is the respondent coefficients (part-worths/utilities)
    Beta = np.zeros((data_dict['L'], data_dict['R']))
    
    for resp in range(data_dict['R']):
        z_resp = 1
        if data_dict['C'] > 1:
            raise NotImplementedError
    
        beta = np.random.multivariate_normal(Gamma, Vbeta)
        if pathology_type:
            beta = pathology(beta, kind=pathology_type)
    
        for scn in range(data_dict['T']):
            X_scn = data_dict['X'][resp, scn]

            U_scn = X_scn.dot(beta) - np.log(-np.log(np.random.uniform(size=data_dict['C'])))
            Y[resp, scn] += np.argmax(U_scn) + 1
    
        Beta[:, resp] += beta

    data_dict['B'] = Beta
    data_dict['Y'] = Y.astype(int)

    return data_dict


def generate_simulated_data(pathology_type="none", use_stan=False):
    if use_stan:
        sm = get_model(model_name='generate_data')
        data_dict = {'A':4, 'L':12, 'T':10, 'R':100, 'C':1}
        data = sm.sampling(data=data_dict,
                           iter=10000,
                           warmup=0,
                           chains=1,
                           refresh=10000,
                           seed=random_seed,
                           algorithm="Fixed_param")
        # pystan fit objects can take a long time to unpack...
        #data_dict.update(data.extract(pars=['X','Y','Z','B']))
        for v in ['X', 'Y', 'Z', 'B']:
            data_dict[v] = data.extract(pars=[v])[v][-1]
        data_dict['B'] = data_dict['B'].T
        data_dict['Y'] = data_dict['Y'].astype(int)
    else:
        data_dict = generate_simulated_design()
        data_dict = compute_beta_response(data_dict, pathology_type=pathology_type)
    return data_dict


def get_model(model_name='mnl_vanilla'):

    with open('../STAN/{0}.stan'.format(model_name), 'r') as f:
        stan_model = f.read()
    
    try:
        sm = pickle.load(open('../STAN/{0}.pkl'.format(model_name), 'rb'))
    
    except:
        sm = pystan.StanModel(model_code=stan_model)
        with open('../STAN/{0}.pkl'.format(model_name), 'wb') as f:
            pickle.dump(sm, f)
    
    return sm


def fit_model_to_data(model, data, sampling_alg=None):
    return model.sampling(
            data,
            iter=niter,
            chains=nchains,
            control={'adapt_delta':.9, 'max_treedepth':treedepth},
            algorithm=sampling_alg,
            init_r=1)



def get_data_dict(pathology_type='none', save=None):

    data_dict = generate_simulated_data(pathology_type=pathology_type)
    
    data_dict['Xtrain'] = data_dict['X'][:nresp_train, :ntask_train, :, :]
    data_dict['Ytrain'] = data_dict['Y'][:nresp_train, :ntask_train]
    data_dict['Xtest'] = data_dict['X'][:nresp_train, -ntask_test:, :, :]
    data_dict['Ytest'] = data_dict['Y'][:nresp_train, -ntask_test:]

    return data_dict


def save_data_dict(data_dict, fpath):
    np.savez(fpath, X=data_dict['X'],
                    Y=data_dict['Y'],
                    Xtrain=data_dict['Xtrain'],
                    Ytrain=data_dict['Ytrain'],
                    Xtest=data_dict['Xtest'],
                    Ytest=data_dict['Ytest'],
                    B=data_dict['B'])


### END ###
