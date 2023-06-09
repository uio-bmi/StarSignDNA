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
from .denovo import denovo as _denovo,cos_sim_matrix



def refit(matrix_file: str, signature_file: str, output_file_exposure: str, output_file_exposure_avg: str, opportunity_file: str = None,
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

    M = pd.read_csv(matrix_file, delimiter='\t').to_numpy().astype(float)
    S = pd.read_csv(signature_file, delimiter=',')
    index_signature=S.index.values.tolist()
    S= S.to_numpy().astype(float)
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
    if data_type == DataType.genome:
        lambd = 0.025
    else:
        lambd = 0.00304
    E, loss, sum_expo = _refit(M, S, O, lambd=lambd)
    E = pd.DataFrame(data=E, columns=index_signature)
    sum_expo = pd.DataFrame(data=sum_expo, columns=index_signature)
    #print(E)
    E.to_csv(output_file_exposure,index=False,header=True,sep='\t')
    sum_expo.to_csv(output_file_exposure_avg,index=False,header=True,sep='\t')


def denovo(matrix_file: str, n_signatures: int, lambd: float,  output_file_exposure: str, output_file_signature: str,
           opportunity_file: str = None, cosmic_file: str= None, max_em_iterations: int = 100,
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
        cos_similarity= cos_sim_matrix(S,cosmic)[0]
        #cos_similarity.to_csv(output_file,index=False,header=True,sep='\t')
        cos_similarity.to_csv("output/cos_sim_skin_pcawg.txt", sep = "\t")
        print(cos_similarity)
    np.savetxt(output_file_exposure, np.array(E))
    np.savetxt(output_file_signature, np.array(S))
    print("--- %s seconds ---" % (time.time() - start_time))


def main():
    app = typer.Typer()
    app.command()(refit)
    app.command()(denovo)
    app()


if __name__ == "__main__":
    main()
