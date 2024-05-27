import numpy as np

from mutational_starsign import denovo
from .test_refit import M, n_signatures, n_samples, n_mutations


def test_denovo_acceptance(M):
    M = M + np.random.rand(n_samples, n_mutations)
    E, S = denovo(M, n_signatures, lambd=0.8)
    assert E.shape == (n_samples, n_signatures)
    assert S.shape == (n_signatures, n_mutations)

    assert E.shape == (n_samples, n_signatures)
