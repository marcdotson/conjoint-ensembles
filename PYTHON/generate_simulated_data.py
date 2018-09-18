import numpy as np

def generate_simulated_design(nresp=100, nscns=10, nalts=4, nlvls=12, ncovs=1):
    
    # beta means
    Gamma = np.random.uniform(-3, 4, size=ncovs * nlvls)
    # beta variance-covariance
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
    
            Y[resp, scn] += np.argmax(U_scn) + 1
            X[resp, scn] += X_scn
    
        Z[:, resp] += z_resp
        Beta[:, resp] += beta.flatten()
    
    # dictionary to store the simulated data and generation parameters
    data_dict = {'X':X,
                 #'Y':Y.astype(int),
                 'Z':Z.T,
                 #'Beta':Beta,
                 'Gamma':Gamma,
                 'Vbeta':Vbeta,
                 'C':nalts,
                 'J':nresp,
                 'G':ncovs,
                 'S':nscns,
                 'K':nlvls}
    return data_dict


def compute_beta(data_dict, pathology=False):

    # Y is the response
    Y = np.zeros((nresp,nscns))
    # Beta is the respondent coefficients (part-worths/utilities)
    Beta = np.zeros((nlvls, nresp))
    
    for resp in range(data_dict['J']):
        z_resp = 1
        if data_dict['G'] > 1:
            raise NotImplementedError
    
        beta = np.random.multivariate_normal(Gamma, Vbeta)
        if pathology:
            raise NotImplementedError
    
        for scn in range(data_dict['S']):
            X_scn = data_dict['X'][resp, scn]
            U_scn = X_scn.dot(beta.flatten()) - np.log(-np.log(np.random.uniform(size=nalts)))
    
            Y[resp, scn] += np.argmax(U_scn) + 1
            X[resp, scn] += X_scn
    
        Z[:, resp] += z_resp
        Beta[:, resp] += beta.flatten()

    data_dict['Beta'] = Beta
    data_dict['Y'] = Y

    pass


def generate_simulated_data():
    data_dict = generate_simulated_design()
    data_dict = compute_beta_response(data_dict)
    return data_dict
