import pandas as pd
import numpy as np

foldername = "01_PathologyNone"

data = np.load("./{0}/data_dict.npz".format(foldername))

R,T,A,L = data['X'].shape
X = data['X'].copy().reshape((R*T*A,L))

Xdf = pd.DataFrame()
Ydf = pd.DataFrame()

X1 = np.repeat([i+1 for i in range(R)], T*A)
X2 = np.repeat([[i+1 for i in range(T)]], R, axis=0).flatten()
X2 = np.repeat(X2, A)
X3 = np.repeat([[i+1 for i in range(A)]], R*T, axis=0).flatten()

Xdf['resp'] = X1
Xdf['task'] = X2
Xdf['alt'] = X3

for l in range(L):
    Xdf['l_{:02d}'.format(l+1)] = X[:,l]


Ydf['resp'] = [i+1 for i in range(R)]
for i in range(T):
    Ydf['task{:02d}'.format(i+1)] = data['Y'][:,i]

utests = np.zeros((100,4))
utests[:,0] = np.random.randint(R)
utests[:,1] = np.random.randint(T)
utests[:,2] = np.random.randint(A)
utests[:,3] = np.random.randint(L)
utests = utests.astype(np.int64)

for i in range(utests.shape[0]):
    r,t,a,l = utests[i]
    df_r = Xdf[Xdf['resp'] == r+1]
    df_rt = df_r[df_r['task'] == t+1]
    df_rta = df_rt[df_rt['alt'] == a+1]
    assert data['X'][r,t,a,l] == df_rta['l_{:02d}'.format(l+1)].iloc[0]
    
    Ydf_r = Ydf[Ydf['resp'] == r+1]
    assert data['Y'][r,t] == Ydf_r['task{:02d}'.format(t+1)].iloc[0]


Xdf.to_csv("./{0}/X.csv".format(foldername), index=False)
Ydf.to_csv("./{0}/Y.csv".format(foldername), index=False)
