#!/usr/bin/python

import sys
from pylmm.input import plink
from pylmm import lmm

from scipy import linalg
import numpy as np
sys.path = ['./src'] + sys.path
from mvLMM import mvLMM


#Y = np.loadtxt("../GAMMA/testData/Y.txt")[:2].T
#K = np.loadtxt("../GAMMA/testData/K.txt")

# geno = np.loadtxt("/common/jungj2/miniconda3/envs/r_env/lib/R/library/gemma2/extdata/mouse100.geno.txt")
Y = np.loadtxt("/common/jungj2/miniconda3/envs/r_env/lib/R/library/gemma2/extdata/mouse100.pheno.txt")
Y = Y[:,[0,5]]
K = np.loadtxt("/common/jungj2/miniconda3/envs/r_env/lib/R/library/gemma2/extdata/mouse100.cXX.txt")
X = np.loadtxt("../rs8275764").reshape(-1,1)


# Create the mvLMM objectq
M = mvLMM(Y, K, norm=True)
# Perform the opitimization
R = M.getMax()

# The important results are then stored in M.mxCor
# Genetic correlation = M.mxCor[0]
# Environmental correlation = M.mxCor[1]
gcor = M.mxCor[0]
ecor = M.mxCor[1]

Psi,Phi = M.getParameterMatrices(gcor, ecor)
print("Psi:")
Psi
print("Phi: ")
Phi
#M.association(X)


# M = mvLMM(Y, K, X0 = X, norm=True)
# R = M.getMax()
# errer!!!
# M.association(X)

###########################################################
### get beta ##############################################
###########################################################

def normPhenos(Y):
    M = Y.shape[1]
    #for i in range(M): Y[:,i] = (Y[:,i] - Y[:,i].mean())/np.sqrt(Y[:,i].var())
    for i in range(M): Y[:,i] = (Y[:,i] - np.mean(Y[:,i]))/np.sqrt(np.var(Y[:,i], ddof=1))
    return Y

def cleanPhenos(Y):
    M = Y.shape[1] 
    for i in range(M):
        y = Y[:,i]
        x = True ^ np.isnan(y)
        if sum(x) == len(y): continue
        m = y[x].mean()
        y[np.isnan(y)] = m
        Y[:,i] = y
    return Y


Kva,Kve = linalg.eigh(K)
leftTransform = Kve.T
with open('leftTransform.npy', 'wb') as f:
    np.save(f, leftTransform)

# Xstar
N = K.shape[0]
X0 = np.ones((N,1))
X0_T = np.dot(leftTransform, X0)
Xstar = np.dot(leftTransform,X)
Xstar = np.hstack([X0_T, np.dot(M.leftTransform,X)])

# Ystar
Y = cleanPhenos(Y)
Y = normPhenos(Y)

Ystar = np.dot(leftTransform, Y)


# Check that they are positive semi-def
Psi_Kva,Psi_Kve = linalg.eigh(Psi)
Phi_Kva,Phi_Kve = linalg.eigh(Phi)

Psi_Kva[Psi_Kva == 0] = 1e-6
Phi_Kva[Phi_Kva == 0] = 1e-6

# Get the D matrix by diagonalizing Psi and Phi 
R = Psi_Kve*np.sqrt(1.0/Psi_Kva) # Now R %*% R.T is = Psi^{-1}
RR = np.dot(np.dot(R.T,Phi),R)
RR_Kva,RR_Kve = linalg.eigh(RR)
D = RR_Kva
rightTransform = np.dot(RR_Kve.T,R.T)


# Get Transformed Y
Yt = np.dot(M.Ystar,rightTransform.T).T.reshape(-1,1)

#X = M.X0_T


# beta_T,mu,beta_T_stderr,_REML_part = M._getBetaT(X,rightTransform,D,Yt)
Ap = rightTransform
P = []
for i in range(M.M): P += (M.Kva + D[i]).tolist()
P = np.array(P)
L = np.kron(Ap, Xstar)


A = L.T * 1.0/(P+1.0)
B = np.dot(A,L)
Bi = linalg.inv(B)
beta_T = np.dot(np.dot(Bi,A),Yt)
beta_T

###########################################################
### get results using new data ############################
###########################################################

Psi = np.array([[2, -0.5],[-0.5, 0.7]])
Phi = np.array([[0.35, 0.38],[0.38, 0.76]])

# Check that they are positive semi-def
Psi_Kva,Psi_Kve = linalg.eigh(Psi)
Phi_Kva,Phi_Kve = linalg.eigh(Phi)

Psi_Kva[Psi_Kva == 0] = 1e-6
Phi_Kva[Phi_Kva == 0] = 1e-6


# Get the D matrix by diagonalizing Psi and Phi 
R = Psi_Kve*np.sqrt(1.0/Psi_Kva) # Now R %*% R.T is = Psi^{-1}
RR = np.dot(np.dot(R.T,Phi),R)
RR_Kva,RR_Kve = linalg.eigh(RR)
D = RR_Kva
rightTransform = np.dot(RR_Kve.T,R.T)


P = []
for i in range(M.M): P += (M.Kva + D[i]).tolist()
P = np.array(P)

# L == Xt
Xt = np.kron(rightTransform, Xstar)
Xt
new_x = (Xt.T * 1.0/(P+1.0)).T
# new_x = np.dot(Xt.T, linalg.inv(np.diag(P))).T

new_y = (Yt.T * 1.0/(P+1.0)).T
# new_y = np.dot(Yt.T, linalg.inv(np.diag(P))).T

t1 = linalg.inv(np.dot(new_x.T, new_x))
t2 = np.dot(new_x.T, new_y)

np.dot(t1, t2)



# from scipy import stats
#res_lm = stats.linregress(X.flatten(), Y[:,0])


