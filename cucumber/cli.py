"""Console script for cucumber."""
import numpy as np
import pandas as pd
from enum import Enum
import time
from numpy import linalg as LA
from scipy import stats
import scipy.spatial as sp
from scipy.optimize import linear_sum_assignment


# todo

class DataType(str, Enum):
    exome = 'exome'
    genome = 'genome'


import typer
from .refit import refit as _refit
from .denovo import denovo as _denovo, cos_sim_matrix
from .count_mutation_cli import count_mutation as count_mutation
from .bootstrapping import bootstrap as _bootstrap


def bootstrap(matrix_file: str, signature_file: str, output_file_exposure_avg: str, output_file_exposure_std: str,
              opportunity_file: str = None, data_type: DataType = DataType.exome):
    M = read_counts(matrix_file)
    S, index_signature = read_signature(signature_file)
    O = read_opportunity(M, opportunity_file)
    lambd = get_lambda(data_type)
    estimated_exposure, exposure_std = _bootstrap(M, S, O, lambd=lambd)
    np.savetxt(output_file_exposure_avg, estimated_exposure, delimiter='\t')
    np.savetxt(output_file_exposure_std, exposure_std, delimiter='\t')


def refit(matrix_file: str, signature_file: str, output_file_exposure: str, output_file_exposure_avg: str,
          opportunity_file: str = None,
          data_type: DataType = DataType.exome):
    '''
    Parameters
    ----------
    matrix_file: str
    signature_file: str
    output_file_exposure: str
    output_file_exposure_avg: str
    opportunity_file: str
    data_type: DataType
    '''

    M = read_counts(matrix_file)
    S, index_signature = read_signature(signature_file)
    O = read_opportunity(M, opportunity_file)
    lambd = get_lambda(data_type)
    E, loss, sum_expo = _refit(M, S, O, lambd=lambd)
    E = pd.DataFrame(data=E, columns=index_signature)
    sum_expo = pd.DataFrame(data=sum_expo, columns=index_signature)
    E.to_csv(output_file_exposure, index=False, header=True, sep='\t')
    sum_expo.to_csv(output_file_exposure_avg, index=False, header=True, sep='\t')


def get_lambda(data_type):
    if data_type == DataType.genome:
        lambd = 0.025
    else:
        lambd = 0.00304
    return lambd


def read_opportunity(M, opportunity_file):
    n_samples = len(M)
    n_mutations = M.shape[1]
    if opportunity_file is not None:
        O = pd.read_csv(opportunity_file, sep='\t', header=None).to_numpy().astype(float)
        O = np.broadcast_to(O, M.shape)
    else:
        O = np.ones((n_samples, n_mutations), dtype=float)
    O = O / np.amax(O).sum(axis=-1, keepdims=True)
    assert O.shape == (n_samples, n_mutations), f'{O.shape} != {(n_samples, n_mutations)}'
    return O


def read_signature(signature_file):
    S = pd.read_csv(signature_file, delimiter=',')
    index_signature = S.index.values.tolist()
    S = S.to_numpy().astype(float)
    return S, index_signature


def read_counts(matrix_file):
    return pd.read_csv(matrix_file, delimiter='\t').to_numpy().astype(float)


def denovo(matrix_file: str, n_signatures: int, lambd: float, output_file_exposure: str, output_file_signature: str,
           opportunity_file: str = None, cosmic_file: str = None, max_em_iterations: int = 100,
           max_gd_iterations: int = 50):
    '''
    Parameters
    ----------
    matrix_file: str
    n_signatures: int
    lambd: float
    output_file_exposure: str
    output_file_signature: str
    opportunity_file: str
    cosmic_file: str
    max_em_iterations
    max_gd_iterations
    '''
    start_time = time.time()
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
        # print(O)
    O = O / np.amax(O).sum(axis=-1, keepdims=True)
    assert O.shape == (n_samples, n_mutations), f'{O.shape} != {(n_samples, n_mutations)}'
    E, S = _denovo(M, n_signatures, lambd, O, em_steps=max_em_iterations, gd_steps=max_gd_iterations)
    if cosmic_file is not None:
        cosmic = pd.read_csv(cosmic_file, delimiter=',')
        cos_similarity = cos_sim_matrix(S, cosmic)[0]
        # cos_similarity.to_csv(output_file,index=False,header=True,sep='\t')
        cos_similarity.to_csv("output/cos_sim_skin_pcawg.txt", sep="\t")
        print(cos_similarity)
    np.savetxt(output_file_exposure, np.array(E))
    np.savetxt(output_file_signature, np.array(S))
    print("--- %s seconds ---" % (time.time() - start_time))


def main():
    app = typer.Typer()
    app.command()(refit)
    app.command()(denovo)
    app.command()(count_mutation)
    app()


if __name__ == "__main__":
    main()
