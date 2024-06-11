import os

import numpy as np
import pytest

from starsigndna.cli import refit, count_mutation
from .test_refit import M, S, O


@pytest.fixture()
def matrix_file(M):
    filename = 'test_data/matrix.csv'
    np.savetxt(filename, (np.array(M)), delimiter=',')
    return filename


@pytest.fixture()
def signature_file(S):
    filename = 'test_data/signatures.csv'
    np.savetxt(filename, (np.array(S)), delimiter=',')
    return filename


@pytest.fixture()
def opportunity_file(O):
    filename = 'test_data/opportunity.csv'
    np.savetxt(filename, (np.array(O)), delimiter=',')
    return filename


@pytest.fixture()
def output_file():
    filename = 'test_data/output.csv'
    if os.path.exists(filename):
        os.remove(filename)
    return filename


def _test_refit_main(matrix_file, signature_file, opportunity_file, output_file):
    refit(matrix_file, signature_file,  opportunity_file)
    #assert os.path.exists(output_file)
    #n_samples = len(M)
    #assert np.loadtxt(output_file).shape == (n_samples, n_signatures)


def test_vcf():
    count_mutation("example_data/few_variants.vcf", "example_data/small_genome.fa", numeric_chromosomes=False)

