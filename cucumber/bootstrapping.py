import numpy as np

from cucumber import refit


def bootstrap_sample(M, rng):
    assert (M.ndim == 2 and M.shape[0] == 1) or M.ndim == 1
    total_count = M.sum()
    ps = M / total_count
    sampled = rng.multinomial(total_count, ps)
    boost = np.broadcast_to(sampled, M.shape)
    # return np.broadcast_to(sampled, M.shape)
    return boost


def bootstrap(M, S, O, n_bootstraps, lambd):
    rng = np.random.default_rng()
    exposures = [refit(bootstrap_sample(M, rng), S, O, lambd=lambd)[0] for _ in range(n_bootstraps)]
    return np.std(exposures, axis = 0)

