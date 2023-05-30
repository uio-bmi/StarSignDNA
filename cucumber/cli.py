"""Console script for cucumber."""
import numpy as np
import pandas as pd
from enum import Enum
# todo

class DataType(str, Enum):
    exome = 'exome'
    genome = 'genome'

import typer
from .refit import refit as _refit


def refit(matrix_file: str, signature_file: str, output_file: str, opportunity_file: str=None, data_type: DataType=DataType.exome):
    '''
    Parameters
    ----------
    matrix_file: str
    signature_file: str
    output_file: str
    opportunity_file: str
    data_type: DataType
    '''

    M = pd.read_csv(matrix_file, delimiter='\t').to_numpy().astype(float)
    S = pd.read_csv(signature_file, delimiter=',').to_numpy().astype(float)
    n_samples = len(M)
    n_signatures = len(S)
    n_mutations = M.shape[1]
    if opportunity_file is not None:
        O = pd.read_csv(opportunity_file, sep='\t', header=None).to_numpy().astype(float)
        O = np.broadcast_to(O, M.shape)
    else:
        O = np.ones((n_samples, n_mutations), dtype=float)
    O = O/np.amax(O).sum(axis=-1, keepdims=True)
    assert O.shape == (n_samples, n_mutations),  f'{O.shape} != {(n_samples, n_mutations)}'
    if data_type == DataType.exome:
        lambd = 0.8
    else:
        lambd = 0.025

    E, loss = _refit(M, S, O, lambd=lambd)
    np.savetxt(output_file, np.array(E))

def denovo(matrix_file, signature_file, n_signatures: int, lambd: float):
    pass


def main():
    app = typer.Typer()
    app.command()(refit)
    app.command()(denovo)
    app()


if __name__ == "__main__":
    main()
