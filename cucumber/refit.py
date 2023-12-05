import numpy as np
import warnings
import math
from scipy.stats import poisson, entropy
from numpy import count_nonzero
from .main_fixed_denovo import Frobinous, Frobinous_reconstuct, running_simulation_refit
from scipy import stats
def refit(M: np.ndarray, S: np.ndarray, O: np.ndarray=None, lambd: float = int, n_iterations: int= 10000) -> np.ndarray:
    '''
    Refit the signatures to the data
    M: Matrix of observed mutational signatures
    S: Mutational Signatures matrix from COSMIC
    O: Matrix of mutational opportunity
    lambd: a penalty parameter that control the rate trade-off between goodness of fit and L-1 norm of the solution
    '''
    # if O is None:
    #     O = np.ones((n_samples, n_mutations), dtype=int)
    # M, S, O = (np.asarray(a) for a in (M, S, O))
    warnings.filterwarnings('ignore')
    np.random.seed(10000)
    n_samples = len(M)
    n_signatures = len(S)
    n_mutations = 96
    if O is None:
        O = np.ones((n_samples, n_mutations), dtype=int)
    M, S, O = (np.asarray(a) for a in (M, S, O))
#    n_mutations = 96
    tmp = np.abs(np.random.laplace(loc=0, scale=1, size=n_samples * n_signatures).reshape(n_samples, n_signatures))
    E = np.full_like(tmp, 0.00001)
    topt = np.float64("Inf")
    # topt= -math.inf
    tedge = np.float64("Inf")
    # tedge = 1
    # print("lambda is", lambd)
    E = running_simulation_refit(E, M, S, O, topt, tedge, lambd, n_steps=n_iterations)
    if (np.any(E < 0)):
        E = np.maximum(E, 0)
    mse_e = Frobinous(M, S, E, O)
    loss = -poisson.logpmf(M, (E @ S) * O)
    # loss[loss == np.float64("Inf")] = 10e+6
    # E[E<=1e-6]=0
    E /= E.sum(axis=-1, keepdims=True)
    E[np.isnan(E)] = 0
    sparsity = 1.0 - ( count_nonzero(E) / float(E.size) )
    E_norm = np.linalg.norm(E, ord=2, axis=None, keepdims=False)
    mse_reconstruct, M_hat = Frobinous_reconstuct(M, S, E, O)
    #print("mse_reconstructed is:", mse_reconstruct)

    print("Sparsity", sparsity)
    print("The MSE is:", np.mean(mse_e))
    # print("PMF of E", np.mean(loss))
    print("E_norm",E_norm)
    print("PMF",np.mean(loss))
    # print("CONV",conv)
    #print("Topt is", topt)
    return E





