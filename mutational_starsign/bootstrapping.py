import numpy as np

from mutational_starsign import refit


def bootstrap_sample_old(M, rng):
    assert (M.ndim == 2 and M.shape[0] == 1) or M.ndim == 1
    total_count = M.sum()
    ps = M / total_count
    sampled = rng.multinomial(total_count, ps)
    boost = np.broadcast_to(sampled, M.shape)
    # return np.broadcast_to(sampled, M.shape)
    return boost


def bootstrap_sample(M):
    #M = M.flatten()
    M = M.astype(int)
    x_new = M.copy()  # Create a copy to avoid modifying the original array
    n = M.sum()
    xc = np.cumsum(M)  # Cumulative sum of x
    # Generate bootstrap sample indices
    k = np.random.choice(range(1, n + 1), size=n, replace=True)  # Use np.random.choice for efficiency
    # Perform bootstrap sampling
    processed = np.zeros(n, dtype=bool)  # Use NumPy array for better performance
    #  for i in range(n):
    for i in range(M.shape[0]):
        ok = np.logical_and(~processed, k <= xc[i])
        x_new[i] = ok.sum()
        processed = np.logical_or(processed, k <= xc[i])
    return x_new


def bootstrap(M,n_bootstraps):
    M = M.flatten()
   # print("MMM", M)
    exposure_boot = np.zeros((n_bootstraps, len(M)))
    bootstrap_M = np.zeros((n_bootstraps, len(M)))
    for b in range(n_bootstraps):
        mboot = bootstrap_sample(M)
      ##  print("AAAA")
        exposure_boot[b,] = mboot  # Store bootstrap samples for exposure (if needed)
        bootstrap_M[b] = mboot
#    print(bootstrap_M.shape)
#    bootstrap_M = float(bootstrap_M)
    bootstrap_M = np.array(bootstrap_M)
 #   exposures = refit(bootstrap_M, S, O, lambd=lambd)[0]
    #   percentile_25 = np.percentile(expo_run, 2.5, axis=0)
    #   percentile_97_5 = np.percentile(expo_run, 97.5, axis=0)
    return bootstrap_M


def bootstrap_old(M, S, O, n_bootstraps, lambd):
    rng = np.random.default_rng()
    expo_run = []
    for _ in range(n_bootstraps):
        exposures = refit(bootstrap_sample(M, rng), S, O, lambd=lambd)[0]
        expo_run.append(exposures)
    percentile_25 = np.percentile(expo_run, 2.5, axis=0)
    percentile_97_5 = np.percentile(expo_run, 97.5, axis=0)
    return percentile_25, percentile_97_5, expo_run
