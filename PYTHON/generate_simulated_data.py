import numpy as np

def pathology(beta, kind=None):
    if kind == 'ANA':
        beta *= np.random.choice([1, 0], size=len(beta), p=[.7, .3])
    elif kind == 'screening':
        beta += np.random.choice([0, -np.inf], size=len(beta), p=[.7, .3])
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
            X_scn = np.random.uniform(size=nalts*nlvls).reshape(nalts,nlvls)
            X[resp, scn] += X_scn
    
        Z[:, resp] += z_resp
    
    # dictionary to store the simulated data and generation parameters
    data_dict = {'X':X,
                 'Z':Z.T,
                 'C':nalts,
                 'J':nresp,
                 'G':ncovs,
                 'S':nscns,
                 'K':nlvls}
    return data_dict


def compute_beta_response(data_dict, pathology_type=None):

    # beta means
    Gamma = np.random.uniform(-3, 4, size=data_dict['G'] * data_dict['K'])
    # beta variance-covariance
    Vbeta = np.diag(np.ones(data_dict['K'])) + .5 * np.ones((data_dict['K'], data_dict['K']))

    # Y is the response
    Y = np.zeros((data_dict['J'], data_dict['S']))
    # Beta is the respondent coefficients (part-worths/utilities)
    Beta = np.zeros((data_dict['K'], data_dict['J']))
    
    for resp in range(data_dict['J']):
        z_resp = 1
        if data_dict['G'] > 1:
            raise NotImplementedError
    
        beta = np.random.multivariate_normal(Gamma, Vbeta)
        if pathology_type:
            beta = pathology(beta, kind=pathology_type)
    
        for scn in range(data_dict['S']):
            X_scn = data_dict['X'][resp, scn]

            U_scn = X_scn.dot(beta.flatten()) - np.log(-np.log(np.random.uniform(size=data_dict['C'])))
            Y[resp, scn] += np.argmax(U_scn) + 1
    
        Beta[:, resp] += beta.flatten()

    data_dict['Beta'] = Beta
    data_dict['Y'] = Y.astype(int)
    data_dict['Gamma'] = Gamma
    data_dict['Vbeta'] = Vbeta
    return data_dict


def generate_simulated_data(pathology_type=None):
    data_dict = generate_simulated_design()
    data_dict = compute_beta_response(data_dict, pathology_type=pathology_type)
    return data_dict
