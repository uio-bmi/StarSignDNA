import numpy as np
import warnings
import math
from scipy.stats import poisson, entropy
from numpy import count_nonzero
from .main_fixed_denovo import Frobinous, Frobinous_reconstuct, running_simulation_refit
from scipy import stats
def refit(M: np.ndarray, S: np.ndarray, O: np.ndarray=None, lambd: float = int, n_iterations: int=1000) -> np.ndarray:
    '''
    Refit the signatures to the data
    M: Matrix of observed mutational signatures
    S: Mutational Signatures matrix from COSMIC
    O: Matrix of mutational opportunity
    lambd: a penalty parameter that control the rate trade-off between goodness of fit and L-1 norm of the solution
    '''
    warnings.filterwarnings('ignore')
    np.random.seed(1000)
    n_samples = len(M)
    n_signatures = len(S)
    n_mutations = 96
    if O is None:
        O = np.ones((n_samples, n_mutations), dtype=int)
    elif O.shape !=M.shape:
        O = np.tile(O, (M.shape[0], 1))
    M, S, O = (np.asarray(a) for a in (M, S, O))
#    n_mutations = 96
    tmp = np.abs(np.random.laplace(loc=0, scale=1, size=n_samples * n_signatures).reshape(n_samples, n_signatures))
    E = np.full_like(tmp, 1) #modify 25 jan
    topt = np.float64("Inf")
    # topt= -math.inf
    tedge = np.float64("Inf")
    expose_mat = []
    pmf_s = []
    mse_s = []
    expo_step = []
    E = running_simulation_refit(E, M, S, O, topt, tedge, lambd, n_steps=n_iterations)
    if (np.any(E < 0)): #comment 25jan
        E = np.maximum(E, 0.00001) #comment 25jan
    mse_e = Frobinous(M, S, E, O)
    loss = -poisson.logpmf(M, (E @ S) * O)
    # print("PMF",np.mean(loss))
    loss[loss == np.float64("Inf")] = 0
    # E[E<=1e-6]=0
    E /= E.sum(axis=-1, keepdims=True)
    E[np.isnan(E)] = 0
    sparsity = 1.0 - ( count_nonzero(E) / float(E.size) )
    E_norm = np.linalg.norm(E, ord='fro', axis=None, keepdims=False)
    mse_reconstruct, M_hat = Frobinous_reconstuct(M, S, E, O)
    print("Sparsity", sparsity)
    print("The MSE is:", np.mean(mse_e))
    # print("PMF of E", np.mean(loss))
    print("E_norm",E_norm)
    print("PMF",np.mean(loss))
    loss_hat = np.mean(loss)
    return E





