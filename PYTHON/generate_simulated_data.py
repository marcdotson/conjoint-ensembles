import pickle
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
    
            Y[resp, scn] += np.argmax(U_scn) + 1
            X[resp, scn] += X_scn
    
        Z[:, resp] += z_resp
        Beta[:, resp] += beta.flatten()
    
    # dictionary to store the simulated data and generation parameters
    data_dict = {'X':X,
                 'Y':Y.astype(int),
                 'Z':Z.T,
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

    with open('./MODELS/HBMNL.stan','r') as f:
        stan_model = f.read()

    try:
        sm = pickle.load(open('./MODELS/HBMNL.pkl', 'rb'))

    except:
        sm = pystan.StanModel(model_code=stan_model)
        with open('./MODELS/HBMNL.pkl','wb') as f:
            pickle.dump(sm, f)

    fit = sm.sampling(data=data_dict, iter=500, chains=4, control={'max_treedepth':5})
    B = fit.extract(pars=['B'])['B'].mean(axis=0)

    plt.figure(figsize=(12,8))
    plt.subplot(311)
    plt.imshow(B.T)
    plt.title("Estimated Betas")
    plt.subplot(312)
    plt.imshow(data_dict['Beta'])
    plt.title("Generated Betas")
    plt.subplot(313)
    plt.plot(np.arange(12), B.T.mean(axis=1), label='Estimated')
    plt.plot(np.arange(12), data_dict['Beta'].mean(axis=1), label='Generated')
    plt.title("Feature Avg Betas")
    plt.show()
