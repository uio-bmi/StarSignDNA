#!/usr/bin/env python3
"""Main module."""
import numpy as np
from scipy.stats import poisson
import main_fixed_denovo
np.random.seed(100)

def main(M, S):
#    M = M.to_numpy()
    S = S.to_numpy()
    n_samples, n_mutations = M.shape
    n_signatures = S.shape[0]
    O = np.ones((n_samples, n_mutations), dtype=int)
    E_true = np.abs(np.random.laplace(loc=0, scale=1, size=n_samples*n_signatures).reshape(n_samples, n_signatures))
    E = np.full_like(E_true, 0.00001)

    topt = np.float64("Inf")
    tedge = np.float64("Inf")
    if(np.any(E < 0)):
        E = np.maximum(E, 0)
    E = main_fixed_denovo.running_simulation_denovo(E, M, S, O, topt, tedge, 0)
    loss = -poisson.logpmf(M, (E@S)*O)
    loss[loss == np.float64("Inf")] = 10e+6
    if(np.any(E < 0)):
        E = np.maximum(E, 0)
    E /= E.sum(axis=-1, keepdims= True)
    E[np.isnan(E)] = 0
    return E, np.mean(loss)
