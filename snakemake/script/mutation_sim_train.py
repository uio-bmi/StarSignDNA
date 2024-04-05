#!/usr/bin/env python3
import numpy as np
from scipy.stats import poisson

import main_fixed_denovo
np.random.seed(200)


def main(M, n_signatures, lambd):
#    M = M.to_numpy()
    n_samples, n_mutations = M.shape
    O = np.ones((n_samples, n_mutations), dtype=int)
    E = np.full((n_samples, n_signatures), 0.00001)
    S = np.random.rand(n_signatures*n_mutations).reshape(n_signatures, n_mutations)

    topt = np.float64("Inf")
    tedge = np.float64("Inf")
    if( np.any(E< 0)):
        E = np.maximum(E,0)
    pmf_e = []
    pmf_s = []
    d_mse_e = []
    d_mse_s = []
    mse_old = np.inf
    pmf_old = np.inf
    for _ in range(25):
        E = main_fixed_denovo.running_simulation_denovo(E, M, S, O, topt, tedge, lambd)# lambd)
        S = main_fixed_denovo.running_simulation_denovo(S.T, M.T, E.T, O.T, topt, tedge, 0).T
        mse_e = main_fixed_denovo.Frobinous(M,S,E,O)
        d_mse_s.append(mse_e)
        loss = -poisson.logpmf(M,(E@S)*O)
        pmf_s.append(np.mean(loss))
        mse_old = mse_e
        pmf_old = np.mean(loss)
    if(np.any(E < 0)):
        E = np.maximum(E, 0)
    E /= E.sum(axis=-1, keepdims= True)
    E[np.isnan(E)] = 0
    if(np.any(S < 0)):
        S = np.maximum(S, 0)
    S /= S.sum(axis=-1, keepdims= True)
    S[np.isnan(S)] = 0
    return np.array(E), np.array(S)
