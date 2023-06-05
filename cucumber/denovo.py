import time

import numpy as np
import pandas as pd
import warnings
from scipy.stats import poisson
from .main_fixed_denovo import Frobinous, running_simulation_new, convergence,cos_sim_matrix
from numpy import linalg as LA
from scipy import stats
from scipy import stats
import scipy.spatial as sp
from scipy.optimize import linear_sum_assignment


def denovo(M: np.ndarray, n_signatures: int, lambd: float, O: np.ndarray = None,em_steps: int = 100,
           gd_steps: int = 50) -> np.ndarray:
    """

    Parameters
    ----------

    M: Matrix of observed mutational signatures
    n_signatures: number of signatures
    lambd: a penalty paramter that control the rate trade-off betweeen goodness of fit and L-1 norm of the solution
    O: Matrix of mutational opportunity
    em_steps:
       Number of steps in the EM algorithm
    gd_steps:
         Number of steps in the gradient descent for E and S
    """

    #warnings.filterwarnings('ignore')
    #time.sleep(10)
    M = np.asanyarray(M)
    n_samples, n_mutations = M.shape
    if O is None:
        O = np.ones((n_samples, n_mutations), dtype=int)
    M, O = (np.asarray(a) for a in (M, O))
    tmp = np.abs(np.random.laplace(loc=0, scale=1, size=n_samples * n_signatures).reshape(n_samples, n_signatures))
    E = np.full_like(tmp, 0.00001)
    S_tmp = np.random.rand(n_signatures, n_mutations).reshape(n_signatures, n_mutations)
    S = np.random.rand(S_tmp.size).reshape(S_tmp.shape)
    topt = np.float64("Inf")
    tedge = np.float64("Inf")
    if np.any(E < 0):
        E = np.maximum(E, 0)
    pmf_e = []
    pmf_s = []
    d_mse_e = []
    d_mse_s = []
    mse_old = np.inf
    pmf_old = np.inf
    conv_iter_1 = 0
    conv_check = 0
    for i in range(em_steps):
        if( np.any(E< 0)):
            E = np.maximum(E,0)
        print("EM step is :",i)
        E = running_simulation_new(E, M, S, O, topt, tedge, lambd, n_steps=gd_steps)  # lambd)
        S = running_simulation_new(S.T, M.T, E.T, O.T, topt, tedge, 0, n_steps=gd_steps).T
        if np.linalg.matrix_rank(S) < S.shape[0]:
            print(np.linalg.matrix_rank(S) < len(S))
            raise Exception("Degenerate Signature Found (rank(S)<len(S)), please provide a lower k")
        mse_e = Frobinous(M, S, E, O)
        d_mse_s.append(mse_e)
        loss = -poisson.logpmf(M, (E @ S) * O)
        pmf_s.append(np.mean(loss))
        mse_old = mse_e
        pmf_old = np.mean(loss)
        print("MSE of S", mse_e)
        print("LOSS S", np.mean(loss))
        conv = convergence(np.mean(loss), pmf_old)
        if conv == True:
            if conv_iter_1 == -1:
                conv_iter_1 = i
                conv_check = 0
            else:
                conv_check = conv_check + 1
        else:
            conv_iter_1 = -1
            conv_check = 0
            if conv_check == 50:
                print("The algorithm converge")
                break
    if (np.any(E < 0)):
        E = np.maximum(E, 0)
    E /= E.sum(axis=-1, keepdims=True)
    E[np.isnan(E)] = 0
    if (np.any(S < 0)):
        S = np.maximum(S, 0)
    S /= S.sum(axis=-1, keepdims=True)
    S[np.isnan(S)] = 0
    return np.array(E), np.array(S)
