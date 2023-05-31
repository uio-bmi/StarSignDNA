import numpy as np
from scipy.stats import poisson, entropy
from .main_fixed_denovo import Frobinous, running_simulation_new, Frobinous_reconstuct
from scipy import stats
def refit(M: np.ndarray, S: np.ndarray, O: np.ndarray, lambd: float = 0.8, n_iterations: int=10) -> np.ndarray:
    '''
    Refit the signatures to the data
    '''
    M, S, O = (np.asarray(a) for a in (M, S, O))
    np.random.seed(10000)
    n_samples = len(M)
    n_signatures = len(S)
    n_mutations = 96
    tmp = np.abs(np.random.laplace(loc=0, scale=1, size=n_samples * n_signatures).reshape(n_samples, n_signatures))
    E = np.full_like(tmp, 0.00001)
    topt = np.float64("Inf")
    tedge = np.float64("Inf")
    if (np.any(E < 0)):
        E = np.maximum(E, 0)
    E = running_simulation_new(E, M, S, O, topt, tedge, lambd, n_steps=n_iterations)
    mse_e = Frobinous(M, S, E, O)
    loss = -poisson.logpmf(M, (E @ S) * O)
    print("MSE of E", mse_e)
    print("PMF of E", np.mean(loss))
    E /= E.sum(axis=-1, keepdims=True)
    E[np.isnan(E)] = 0
    sum_expo = E.sum(axis=0, keepdims=True) / len(E)

    mse_reconstruct, M_hat = Frobinous_reconstuct(M, S, E, O)
    print("mse_reconstructed is:", mse_reconstruct)
    entropy_reco = entropy(M, axis=1, qk=M_hat)
    entropy_reco_mean = np.mean(entropy_reco)
    print("The relative entropy is:", entropy_reco_mean)
    rho, pval = stats.spearmanr(np.transpose(M), np.transpose(M_hat))
    print("The speraman coefficient is:", np.mean(rho), np.mean(pval))
    print("The final log poisson PMF is:", np.mean(loss))
    return E, np.mean(loss)




