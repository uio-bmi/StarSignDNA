import numpy as np
from scipy.stats import poisson
from .main_fixed_denovo import Frobinous, running_simulation_new, convergence


def denovo(M: np.ndarray, n_signatures: int, lambd: float, em_steps: int=10, gd_steps: int =50):
    """

    Parameters
    ----------
    M
    n_signatures
    lambd
    em_steps:
       Number of steps in the EM alogrithm
    gd_steps:
         Number of steps in the gradient descent for E and S
    """
    n_samples, n_mutations = M.shape
    # O = np.ones((n_samples, n_mutations), dtype=int)
    E = np.full((n_samples, n_signatures), 0.00001)
    S = np.random.rand(n_signatures * n_mutations).reshape(n_signatures, n_mutations)
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
    for _ in range(em_steps):
        E = running_simulation_new(E, M, S, O, topt, tedge, lambd, n_steps=gd_steps)  # lambd)
        S = running_simulation_new(S.T, M.T, E.T, O.T, topt, tedge, 0, n_steps=gd_steps).T
        mse_e = Frobinous(M, S, E, O)
        d_mse_s.append(mse_e)
        loss = -poisson.logpmf(M, (E @ S) * O)
        pmf_s.append(np.mean(loss))
        mse_old = mse_e
        pmf_old = np.mean(loss)
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
