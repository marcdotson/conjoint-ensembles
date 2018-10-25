'''
Functions in this file:

    pathology(beta, kind=None)

    generate_simulated_design(nresp, nscns, nalts, nlvls, ncovs)

    compute_beta_response(data_dict, pathology_type=None)

    generate_simulated_data(pathology_type=None)

    fit_model(data_dict, model_name)

    get_loo_list(fit)

'''


from . import psis
import pystan
import pickle
import numpy as np
from scipy.optimize import minimize, LinearConstraint


def pathology(beta, kind=None, prob=[.5, .5]):
    if kind == 'ANA':
        beta *= np.random.choice([1, 0], size=len(beta), p=prob)
    elif kind == 'screening':
        beta += np.random.choice([0, -np.inf], size=len(beta), p=prob)
    elif kind == 'systematicANA':
        pathology_vector = np.ones_like(beta)
        pathology_vector[:int(len(beta)//2)] = 0
        beta *= pathology_vector
    elif kind == 'randomANA':
        if int(np.random.choice([1, 0], p=prob)):
            pathology_vector = np.ones_like(beta)
            pathology_vector[:int(len(beta)//2)] = 0
            beta *= pathology_vector
    return beta


def generate_simulated_design(nresp=100, nscns=10, nalts=4, nlvls=12, ncovs=1):
    
    # X is the experimental design
    X = np.zeros((nresp, nscns, nalts, nlvls))
    # Z is the covariates
    Z = np.zeros((ncovs, nresp))
    
    for resp in range(nresp):
        z_resp = 1
        if ncovs > 1:
            raise NotImplementedError
    
        for scn in range(nscns):
            X_scn = np.random.choice([0,1], p =[.5,.5], size=nalts*nlvls).reshape(nalts,nlvls)
            X[resp, scn] += X_scn
    
        Z[:, resp] += z_resp
    
    # dictionary to store the simulated data and generation parameters
    data_dict = {'X':X,
                 'Z':Z.T,
                 'A':nalts,
                 'R':nresp,
                 'C':ncovs,
                 'T':nscns,
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

            U_scn = X_scn.dot(beta.flatten()) - np.log(-np.log(np.random.uniform(size=data_dict['C'])))
            Y[resp, scn] += np.argmax(U_scn) + 1
    
        Beta[:, resp] += beta.flatten()

    data_dict['Beta'] = Beta
    data_dict['Y'] = Y.astype(int)
    data_dict['Gamma'] = Gamma
    data_dict['Vbeta'] = Vbeta

    if pathology_type == 'ANA':
        data_dict['w'] = np.random.binomial(1, .5, size=data_dict['T'])

    return data_dict


def generate_simulated_data(pathology_type=None):
    data_dict = generate_simulated_design()
    data_dict = compute_beta_response(data_dict, pathology_type=pathology_type)
    return data_dict




def fit_model(data_dict, model_name='HBMNL_vanilla'):

    with open('./STAN/{0}.stan'.format(model_name), 'r') as f:
        stan_model = f.read()
    
    try:
        sm = pickle.load(open('./STAN/{0}.pkl'.format(model_name), 'rb'))
    
    except:
        sm = pystan.StanModel(model_code=stan_model)
        with open('./STAN/{0}.pkl'.format(model_name), 'wb') as f:
            pickle.dump(sm, f)
    
    fit = sm.sampling(data=data_dict, iter=1000, chains=2)

    return fit

#def get_loos(fit):
#
#    log_lik = fit.extract(pars='log_lik')['log_lik']
#
#    LL = np.zeros((log_lik.shape[0], log_lik.shape[1]*log_lik.shape[2]))
#    for i in range(log_lik.shape[0]):
#        LL[i] = log_lik[i].flatten()
#
#    return psis.psisloo(LL)
#
#
#def stacking_weights(LL_list):
#    lpd_point = np.vstack(LL_list).T
#    N = lpd_point.shape[0]
#    K = lpd_point.shape[1]
#    exp_lpd_point = np.exp(lpd_point)
#
#    # neg_log_score_loo
#    def f(w):
#        w_full = np.hstack((w, 1-np.sum(w)))
#        S = 0
#        for i in range(N):
#            S += np.log(np.exp(lpd_point[i, :]).dot(w_full))
#        return -S
#
#    # grad_neg_log_score_loo
#    def grad_f(w):
#        w_full = np.hstack((w, 1-np.sum(w)))
#        grad = np.zeros(K-1)
#        for k in range(K-1):
#            for i in range(N):
#                grad[k] += (exp_lpd_point[i,k] - exp_lpd_point[i,-1]) / (exp_lpd_point[i, :].dot(w_full))
#        return -grad
#
#    ui = np.vstack((-np.ones(K-1), np.diag(np.ones(K-1))))
#    ci = np.zeros(K)
#    ci[0] = -1
#    x0 = (1/K)*np.ones(K-1)
#    lincon = ({'type': 
#    out = minimize(f, x0, method='COBYLA', jac=grad_f, constraints=[lincon])
#    return out
#
