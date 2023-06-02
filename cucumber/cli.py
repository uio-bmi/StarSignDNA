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
from .denovo import denovo as _denovo


def refit(matrix_file: str, signature_file: str, output_file: str, opportunity_file: str = None,
          data_type: DataType = DataType.exome):
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
    O = O / np.amax(O).sum(axis=-1, keepdims=True)
    assert O.shape == (n_samples, n_mutations), f'{O.shape} != {(n_samples, n_mutations)}'
    if data_type == DataType.exome:
        lambd = 0.8
    else:
        lambd = 0.025

    E, loss = _refit(M, S, O, lambd=lambd)
    np.savetxt(output_file, np.array(E))


def denovo(matrix_file: str, n_signatures: int, lambd: float, output_file_exposure: str, output_file_signature: str,
           opportunity_file: str = None, max_em_iterations: int = 10, max_gd_iterations: int = 50):
    '''
    Parameters
    ----------
    matrix_file: str
    n_signatures: int
    lambd: float
    output_file_exposure: str
    output_file_signature: str
    opportunity_file: str
    max_em_iterations
    max_gd_iterations
    '''

    M = pd.read_csv(matrix_file, delimiter='\t').to_numpy().astype(float)
    n_samples = len(M)
    n_signatures = n_signatures
    lambd = lambd
    n_mutations = M.shape[1]
    if opportunity_file is not None:
        O = pd.read_csv(opportunity_file, sep='\t', header=None).to_numpy().astype(float)
        O = np.broadcast_to(O, M.shape)
    else:
        O = np.ones((n_samples, n_mutations), dtype=float)
        print(O)
    O = O / np.amax(O).sum(axis=-1, keepdims=True)
    assert O.shape == (n_samples, n_mutations), f'{O.shape} != {(n_samples, n_mutations)}'
    E, S = _denovo(M, O, n_signatures, lambd, em_steps=max_em_iterations, gd_steps=max_gd_iterations)
    np.savetxt(output_file_exposure, np.array(E))
    np.savetxt(output_file_signature, np.array(S))


def main():
    app = typer.Typer()
    app.command()(refit)
    app.command()(denovo)
    app()


if __name__ == "__main__":
    main()
