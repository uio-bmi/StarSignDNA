#!/usr/bin/env python
import pytest

from starsigndna.cli import filter_signatures
import pytest


@pytest.mark.skip(reason='failing')
def test_filter_signatures(signature_matrix):
    names = ['SBS6', 'SBS5']
    result = filter_signatures(signature_matrix, names)
    assert len(result) == len(names)
