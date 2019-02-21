import pystan
import pickle
import numpy as np
import pandas as pd


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


def generate_simulated_data(pathology_type="none"):
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


def fit_model_to_data(model, data, **kwargs):
    """Runs the Stan sampler for model on data according to kwargs."""
    return model.sampling(data, **kwargs)


def get_data(path_to_data, holdout=5):
    """
    1. get X.csv and Y.csv
    2. reformat data into ndarrays
    3. split into Xtrain, Ytrain, Xtest, Ytest (take out some choice tasks)
    4. save as dictionary of arrays X,Y,Xtrain,Ytrain,Xtest,Ytest
    """

    # dictionary to store the ndarrays
    data_dict = dict()

    # read in csv files and reformat
    Xdf = pd.read_csv(path_to_data + "X.csv")
    Ydf = pd.read_csv(path_to_data + "Y.csv")

    # determine data dimensions
    nresp = Xdf['resp'].max()
    ntask = Xdf['task'].max()
    nalts = Xdf['alt'].max()
    nlvls = len(Xdf.columns) - 3

    # format X and Y
    X = Xdf.drop(['resp','task','alt'],axis=1).values
    Y = Ydf.drop(['resp'],axis=1).values
    X = X.reshape((nresp,ntask,nalts,nlvls))
    Y = Y.astype(np.int64).reshape((nresp,ntask))

    # store X and Y in dictionary
    data_dict['X'] = X
    data_dict['Y'] = Y

    # store training and test sets
    data_dict['Xtrain'] = X[:,:-holdout,:,:]
    data_dict['Ytrain'] = Y[:,:-holdout]
    data_dict['Xtest'] = X[:,-holdout:,:,:]
    data_dict['Ytest'] = Y[:,-holdout:]

    return data_dict



### END ###
