import numpy as np

from cucumber import refit


def bootstrap_sample(M):
    assert (M.ndim == 2 and  M.shape[0] == 1)  or M.ndim ==1
    ps = M/M.sum
    rng = np.random.default_rng()
    sampled = rng.multinomial(M.sum(), ps)
    return np.broadcast_to(sampled, M.shape)


def bootstrap(M, S, O,  n_bootstraps = 100, lambd=0):
    exposures = [refit(bootstrap_sample(M), S, O, lambd=lambd) for _ in range(n_bootstraps)]
    return np.mean(exposures, axis=0), np.std(exposures, axis=0)
