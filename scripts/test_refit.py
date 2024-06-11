import numpy as np
import pytest

from starsigndna import denovo

n_samples = 2
n_signatures = 3
n_mutations = 4


@pytest.fixture
def M():
    return np.ones((n_samples, n_mutations))


@pytest.fixture
def S():
    return np.ones((n_signatures, n_mutations))


@pytest.fixture
def O():
    return np.ones((n_samples, n_mutations))


def test_refit_acceptance(M):
    E, S = denovo(M, n_signatures, lambd=0.8)
    assert E.shape == (n_samples, n_signatures)
    assert M.shape == (n_samples, n_mutations)

