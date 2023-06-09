import numpy as np

from cucumber import refit


def bootstrap_sample(M, rng):
    assert (M.ndim == 2 and  M.shape[0] == 1)  or M.ndim ==1
    total_count = M.sum()
    ps = M / total_count
    sampled = rng.multinomial(total_count, ps)
    return np.broadcast_to(sampled, M.shape)


def bootstrap(M, S, O,  n_bootstraps = 100, lambd=0):
    rng = np.random.default_rng()
    exposures = [refit(bootstrap_sample(M, rng), S, O, lambd=lambd) for _ in range(n_bootstraps)]
    return np.mean(exposures, axis=0), np.std(exposures, axis=0)
