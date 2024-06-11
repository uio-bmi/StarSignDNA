from pathlib import Path

import pytest
from starsigndna.cli import read_signature



@pytest.fixture()
def data_path():
    return Path(__file__).parent.parent / 'example_data'

@pytest.fixture
def signature_matrix(data_path):
    return read_signature(data_path / 'COSMICv34.txt')
