import pandas as pd
import numpy as np

fpath = "./DATA/03_ArtificialFlowers/"
design_fname = "coded.design.csv"
survey_fname = "ecoflower.data.csv"
T = 6

conjoint_columns = ["0.Version"]
for i in range(T):
    conjoint_columns.append("C{0}".format(i+1))



data_df = pd.read_csv(fpath + survey_fname) #, encoding='ISO-8859-1')
Ydf = data_df[conjoint_columns]
Xdf = pd.read_csv(fpath + design_fname)
#Xdf = Xdf[Xdf.columns[1:]]

print(Ydf.head())
print(Xdf.head())
input("continue?")

R = Ydf.shape[0]
T = Ydf.shape[1] - 1
A = 4
L = Xdf.shape[1] - 3
C = 1

X = np.zeros((R, T, A, L))
Y = np.zeros((R,T))
print(X.shape)

for r in range(R):
    version = Ydf.iloc[r]['0.Version']
    Xdesign = Xdf[Xdf['Version'] == version]
    for t in range(T):
        Xtask = Xdesign[Xdesign['Task'] == t+1]
        A_L = Xtask.drop(['Version','Task','Concept'],axis=1).values
        X[r,t] = A_L
        Y[r,t] = Ydf.iloc[r]['C{0}'.format(t+1)]

print(Y)
print(X.shape)
print(Y.shape)

input("PRESS ENTER TO SAVE ({0})".format(fpath+"data_dict.npz"))

np.savez(fpath + "data_dict.npz",
         X=X,
         Y=Y,
         Xtrain=X[:,:4,:,:],
         Ytrain=Y[:,:4],
         Xtest=X[:,4:,:,:],
         Ytest=Y[:,4:])

