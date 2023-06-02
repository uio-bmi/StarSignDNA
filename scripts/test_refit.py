import numpy as np
import pytest

from cucumber.refit import refit

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


def test_refit_acceptance(M, S, O):
    E, loss = refit(M, S, O)
    assert E.shape == (n_samples, n_signatures)
