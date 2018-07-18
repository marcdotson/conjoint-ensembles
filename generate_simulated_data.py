import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pystan

def generate_simulated_data():
    
    nresp = 100
    nscns = 10
    nalts = 4
    nlvls = 12
    ncovs = 1
    
    Gamma = np.random.uniform(-3, 4, size=ncovs * nlvls)
    Vbeta = np.diag(np.ones(nlvls)) + .5 * np.ones((nlvls, nlvls))
    
    # Y is the response
    Y = np.zeros((nresp,nscns))
    # X is the experimental design
    X = np.zeros((nresp, nscns, nalts, nlvls))
    # Z is the covariates
    Z = np.zeros((ncovs, nresp))
    # Beta is the respondent coefficients (part-worths/utilities)
    Beta = np.zeros((nlvls, nresp))
    
    for resp in range(nresp):
        z_resp = 1
        if ncovs > 1:
            raise NotImplementedError
    
        beta = np.random.multivariate_normal(Gamma, Vbeta)
    
        for scn in range(nscns):
            X_scn = np.random.uniform(size=nalts*nlvls).reshape(nalts,nlvls)
            U_scn = X_scn.dot(beta.flatten()) - np.log(-np.log(np.random.uniform(size=nalts)))
    
            Y[resp, scn] += np.argmax(U_scn)
            X[resp, scn] += X_scn
    
        Z[:, resp] += z_resp
        Beta[:, resp] += beta.flatten()
    
    # dictionary to store the simulated data and generation parameters
    data_dict = {'X':X,
                 'Y':Y,
                 'Z':Z,
                 'Beta':Beta,
                 'Gamma':Gamma,
                 'Vbeta':Vbeta,
                 'C':nalts,
                 'J':nresp,
                 'G':ncovs,
                 'S':nscns,
                 'K':nlvls}
    return data_dict

if __name__ == '__main__':
    data_dict = generate_simulated_data()

    sm = pystan.StanModel(model_code='./HBMNL.stan')
    fit = sm.sampling(data=data_dict, iter=300, chains=2, max_tree_depth=3)
