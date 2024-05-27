import pytest

from starsign.cli import filter_signatures


def test_filter_signatures(signature_matrix):
    names = ['SBS6', 'SBS5']
    result = filter_signatures(signature_matrix, names)
    assert len(result) == len(names)
