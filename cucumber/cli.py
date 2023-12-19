"""Console script for cucumber."""
import numpy as np
import pandas as pd
from enum import Enum
import time
import os
from matplotlib import pyplot as plt
import seaborn as sns
import string
import multiprocessing

np.random.seed(1000)


# todo


class DataType(str, Enum):
    exome = 'exome'
    genome = 'genome'


import typer
from .refit import refit as _refit
from .denovo import denovo as _denovo, cos_sim_matrix
from .count_mutation_cli import count_mutation
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


# def plot(file, output_folder='output/'):
#     plt.style.use('default')
#     file.plot(kind="bar", figsize=(10, 4))
#     plt.xticks(fontsize=8)
#     plt.yticks(fontsize=8)
#     plt.ylabel("Mutation fraction")
#     plt.xlabel("Signatures")
#     plt.tight_layout()
#     return plt.savefig(f"{output_folder}/plot_cohort.png", dpi=600)


def single_plot(file):
    file = file.transpose()
    file.columns = ['Signatures', 'std_dev']
    # filtered_data = file[file['Signatures'] != 0]
    filtered_data = file[file['Signatures'] >= 0.05]

    # Set up the figure and axis
    plt.style.use('default')
    fig, ax = plt.subplots()
    fig, ax = plt.subplots(figsize=(20, 8))

    # Plot the bar plot with error bars for non-zero values
    ax.bar(filtered_data.index, filtered_data['Signatures'], yerr=filtered_data['std_dev'],
           color='blue', alpha=0.5, align='center', capsize=3)
    plt.ylim(0, None)
    plt.xticks(rotation=45)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.xticks(rotation=45)
    # Set labels and title
    ax.set_xlabel('Signature exposures')
    ax.set_ylabel('Mutation fraction')
    ax.set_title('Single Sample Mutational Signatures')
    return plt


def cohort_plot(file):
    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(20, 8))
    file.plot(kind="bar", figsize=(10, 4))
    plt.xticks(rotation=45)
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=7)
    plt.ylabel("Mutation fraction")
    plt.xlabel("Signature exposures")
    plt.tight_layout()
    return plt


def cohort_violin(file):
    fig, ax = plt.subplots(figsize=(20, 8))
    sns.violinplot(data=file, ax=ax, scale="count", cut=0)
    plt.xticks(rotation=45)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.ylabel("Mutation fraction")
    plt.xlabel("Signature exposures")
    return plt


def plot_profile(data):
    plt.style.use('default')
    # data.drop(columns=data.columns[0], axis=1, inplace=True)
    data = data.T
    header = data.index
    data_list = data.values.tolist()
    len(data_list)
    mutation_categories = np.arange(1, 97)
    mutation_labels = [
        'C>A', 'C>A', 'C>A', 'C>A', 'C>A', 'C>A', 'C>A', 'C>A',
        'C>A', 'C>A', 'C>A', 'C>A', 'C>A', 'C>A', 'C>A', 'C>A',
        'C>G', 'C>G', 'C>G', 'C>G', 'C>G', 'C>G', 'C>G', 'C>G',
        'C>G', 'C>G', 'C>G', 'C>G', 'C>G', 'C>G', 'C>G', 'C>G',
        'C>T', 'C>T', 'C>T', 'C>T', 'C>T', 'C>T', 'C>T', 'C>T',
        'C>T', 'C>T', 'C>T', 'C>T', 'C>T', 'C>T', 'C>T', 'C>T',
        'T>A', 'T>A', 'T>A', 'T>A', 'T>A', 'T>A', 'T>A', 'T>A',
        'T>A', 'T>A', 'T>A', 'T>A', 'T>A', 'T>A', 'T>A', 'T>A',
        'T>C', 'T>C', 'T>C', 'T>C', 'T>C', 'T>C', 'T>C', 'T>C',
        'T>C', 'T>C', 'T>C', 'T>C', 'T>C', 'T>C', 'T>C', 'T>C',
        'T>G', 'T>G', 'T>G', 'T>G', 'T>G', 'T>G', 'T>G', 'T>G',
        'T>G', 'T>G', 'T>G', 'T>G', 'T>G', 'T>G', 'T>G', 'T>G'
    ]

    color_groups = {
        'C>A': 'red',
        'C>G': 'blue',
        'C>T': 'green',
        'T>A': 'orange',
        'T>C': 'purple',
        'T>G': 'brown'
    }

    mutation_colors = [color_groups[label] for label in mutation_labels]

    mutation_matrix = np.zeros((1, 96))
    fig = plt.figure(figsize=(15, 8))
    n_rows = len(data_list)
    fig_rows = round(n_rows / 2, 0)
    # print("NNN", fig_rows)
    fi_rows_rest = len(data_list) % 2
    # print("RRR", fi_rows_rest)
    for id_x in range(len(data_list)):
        plt.subplot(fig_rows + fi_rows_rest, 2, id_x + 1)
        plt.subplots_adjust(left=0.1,
                            bottom=0.1,
                            right=0.9,
                            top=0.9,
                            wspace=0.4,
                            hspace=0.4)
        plt.bar(mutation_categories, data_list[id_x], color=mutation_colors)
        plt.xlabel('Mutation Categories')
        plt.ylabel('Mutations')
        plt.title(f'{header[id_x]}')
        plt.xticks(mutation_categories[::16], mutation_labels[::16], ha='center')
        plt.xlim(1, 97)
        plt.grid(b=True, color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)
        legend_handles = [plt.Rectangle((0, 0), 1, 1, color=color) for color in color_groups.values()]
        legend_labels = color_groups.keys()
        plt.legend(legend_handles, legend_labels)
    return plt


def refit(matrix_file: str, signature_file: str,
          opportunity_file: str = None, ref_genome: str = None,
          data_type: DataType = DataType.exome, n_bootstraps: int = 25, numeric_chromosomes: bool = False,
          genotyped: bool = True, output_folder: str = 'output/', cancer_type: str = None):
    '''
    Parameters
    ----------
    numeric_chromosomes
    n_bootstraps
    matrix_file: str
    signature_file: str
    output_file_exposure: str
    opportunity_file: str
    data_type: DataType
    numeric_chromosomes: bool
        True if chromosome names in vcf are '1', '2', '3'. False if 'chr1', 'chr2', 'chr3'
    genotyped: bool
        True if the VCF file has genotype information for many samples
    '''
    start_time = time.time()
    # reference genome
    # ref_genome = '/Users/bope/Documents/MutSig/scientafellow/packages/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
    file_name, file_extension = os.path.splitext(matrix_file)
    if file_extension == '.vcf':
        assert ref_genome is not None, 'Please provide a reference genome along with the vcf file'
        count_mutation(matrix_file, ref_genome, f'{output_folder}/matrix.csv', numeric_chromosomes, genotyped)
        matrix_file = f'{output_folder}/matrix.csv'
    #M = read_counts(matrix_file)
    M, index_matrix = read_counts(matrix_file)
    #print(index_matrix)
    #M = M[0:5,0:96]
    #    S, index_signature = read_signature(signature_file)
    S, index_signature = read_signature(signature_file)
    if cancer_type is not None:
        if cancer_type == 'bcla':
            index = [0,1,3,4 ,18]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'brca':
            index = [0,1,2,4,18,24,37,41,50]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'chol':
            index = [0,1,4,18,24,47,48,49,67]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'gbm':
            index = [0,4,37,47,48,49]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'lgg':
            index = [0,4]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'cesc':
            index = [0,1,4,18]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'coad':
            index = [0,4,20,47,48,49,53]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'esca':
            index = [0,2,4,18,22,23,24, 47,48,49]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'uvm':
            index = [0,4,47,48,49,60,61]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'hnsc':
            index = [0,1,3,4,18,47,48,49,54]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'kich':
            index = [0,1,18,36, 47,48,49]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'kirp':
            index = [0,1,4,18,54,58]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'kirc':
            index = [0,4,28,29,47,48,49,50]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'lihc':
            index = [0,2,3,4,17,21,24,28,29, 36,47,48,49]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'luad':
            index = [0,1,3,4,18,47,48,54]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'lusc':
            index = [0,1,3,4,18,54]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'dlbc':
            index = [1,2,4,5,11,18,22,23,41,43,44,47,48,49,65]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'dlbc':
            index = [0,1,2,4,5, 11,18,22,23,41,43,44,47,48,49]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'laml':
            index = [0,4,52]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'ov':
            index = [0,1,2,4,18,24,47,48,49]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'paad':
            index = [0,1,2,4,18,22,23,24,47,48,49]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'prad':
            index = [0,2,4,24,47,48,49,67]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'sarc':
            index = [0,4,6,7,8,9,24,68]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'skcm':
            index = [0,4,6,7,8,9,45,47,48,49,54,58]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'stad':
            index = [0,1,2,18,20,22,23,24,26,47,48,49]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'thca':
            index = [0,1,4,18,47,48,49,52,54,58,62]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'ucec':
            index = [0,1,4,12,13,14,15,18,19,20,35,53]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        # elif cancer_type == 'test':
        #     index = [0,1,2,3,4,5,6,7,8,9,10]
        #     S = S[index]
        #     index_signature = [index_signature[i] for i in index]
        else:
            raise ValueError(f'Unknown cancer type {cancer_type}')
    O = read_opportunity(M, opportunity_file)
    #lambd = get_lambda(data_type)
    lambd = 5 * np.sqrt(len(M)/100)
    if (M.ndim == 2 and M.shape[0] == 1) or M.ndim == 1:
        ##normaize_m
        # M /= M.sum(axis=-1, keepdims=True)
        #print(lambd)
        E = _refit(M, S, O, lambd=lambd)[0]
        E_std = _bootstrap(M, S, O, n_bootstraps, lambd=lambd)
        E = [E, E_std]
        E = pd.DataFrame(data=E, columns=index_signature, index=['Signature', 'std_dev'])
        #   E.drop(columns=E.columns[0], axis=1, inplace=True)
        plot = single_plot(E)
        # E.to_csv(output_file_exposure, index=True, header=True, sep='\t')
        E.to_csv(f"{output_folder}/exposure_signature.txt", index=True, header=True, sep='\t')
        plot.savefig(f"{output_folder}/plot_exposure_signature.png", dpi=600)

    else:
        #print(lambd)
        E = _refit(M, S, O, lambd=lambd)
        #print(O)
        sum_expo = E.sum(axis=0, keepdims=True) / len(E)
        E = pd.DataFrame(data=E, columns=index_signature, index=index_matrix)
 #       E.to_csv(f'{output_folder}/cucumber_cohort_exposure.txt', index=False, header=True, sep='\t')
        E.to_csv(f'{output_folder}/cucumber_cohort_exposure.txt', index=index_matrix, header=True, sep='\t')
        sum_expo = pd.DataFrame(data=sum_expo, columns=index_signature, index=['Signatures'])
        # Add the row index to the DataFrame
        #sum_expo = sum_expo.index_matrix()
        # print(sum_expo)
        sum_expo = np.transpose(sum_expo)
        plot_summary = cohort_plot(sum_expo)
        plot_summary.savefig(f"{output_folder}/exposures_cohort_dotplot.png", dpi=600)
        sort_E = sum_expo.sort_values(by=['Signatures'], ascending=False)
        sort_E = sort_E.iloc[:5, 0:]
        plot_top_five = cohort_plot(sort_E)
        plot_top_five.savefig(f"{output_folder}/exposures_cohort_top_5.png", dpi=600)
        plot_variance = cohort_violin(E)
        plot_variance.savefig(f"{output_folder}/exposures_cohort_variance.png", dpi=600)
        # sum_expo.to_csv(output_file_exposure_avg, index=index_signature, header=T67.png", dpi=600)
        # sum_expo.to_csv(output_file_exposure_avg, index=index_signature, header=True, sep='\t')
        sum_expo.to_csv(f'{output_folder}/average_exposure.txt', index=index_signature, header=True,
                        sep='\t')
    print("--- %s seconds ---" % (time.time() - start_time))


# def get_lambda(data_type):
#     if data_type == DataType.genome:
#         lambd = 10
#     else:
#         lambd = 5
#     return lambd


def read_opportunity(M, opportunity_file):
    n_samples = len(M)
    n_mutations = M.shape[1]
    if opportunity_file is not None:
        if opportunity_file == "human-genome":
            O = np.array([1.14e+08, 6.60e+07, 1.43e+07, 9.12e+07, # C>A @ AC[ACGT]
                  1.05e+08, 7.46e+07, 1.57e+07, 1.01e+08, # C>A @ CC[ACGT]
                  8.17e+07, 6.76e+07, 1.35e+07, 7.93e+07, # C>A @ GC[ACGT]
                  1.11e+08, 8.75e+07, 1.25e+07, 1.25e+08, # C>A @ TC[ACGT]
                  1.14e+08, 6.60e+07, 1.43e+07, 9.12e+07, # C>G @ AC[ACGT]
                  1.05e+08, 7.46e+07, 1.57e+07, 1.01e+08, # C>G @ CC[ACGT]
                  8.17e+07, 6.76e+07, 1.35e+07, 7.93e+07, # C>G @ GC[ACGT]
                  1.11e+08, 8.75e+07, 1.25e+07, 1.25e+08, # C>G @ TC[ACGT]
                  1.14e+08, 6.60e+07, 1.43e+07, 9.12e+07, # C>T @ AC[ACGT]
                  1.05e+08, 7.46e+07, 1.57e+07, 1.01e+08, # C>T @ CC[ACGT]
                  8.17e+07, 6.76e+07, 1.35e+07, 7.93e+07, # C>T @ GC[ACGT]
                  1.11e+08, 8.75e+07, 1.25e+07, 1.25e+08, # C>T @ TC[ACGT]
                  1.17e+08, 7.57e+07, 1.04e+08, 1.41e+08, # T>A @ AC[ACGT]
                  7.31e+07, 9.55e+07, 1.15e+08, 1.13e+08, # T>A @ CC[ACGT]
                  6.43e+07, 5.36e+07, 8.52e+07, 8.27e+07, # T>A @ GC[ACGT]
                  1.18e+08, 1.12e+08, 1.07e+08, 2.18e+08, # T>A @ TC[ACGT]
                  1.17e+08, 7.57e+07, 1.04e+08, 1.41e+08, # T>C @ AC[ACGT]
                  7.31e+07, 9.55e+07, 1.15e+08, 1.13e+08, # T>C @ CC[ACGT]
                  6.43e+07, 5.36e+07, 8.52e+07, 8.27e+07, # T>C @ GC[ACGT]
                  1.18e+08, 1.12e+08, 1.07e+08, 2.18e+08, # T>C @ TC[ACGT]
                  1.17e+08, 7.57e+07, 1.04e+08, 1.41e+08, # T>G @ AC[ACGT]
                  7.31e+07, 9.55e+07, 1.15e+08, 1.13e+08, # T>G @ AC[ACGT]
                  6.43e+07, 5.36e+07, 8.52e+07, 8.27e+07, # T>G @ AG[ACGT]
                  1.18e+08, 1.12e+08, 1.07e+08, 2.18e+08]) # T>G @ AT[ACGT]])
        elif opportunity_file == "human-exome":
            O = np.array([1940794, 1442408, 514826, 1403756,
                  2277398, 2318284, 774498, 2269674,
                  1740752, 1968596, 631872, 1734468,
                  1799540, 1910984, 398440, 2024770,
                  1940794, 1442408, 514826, 1403756,
                  2277398, 2318284, 774498, 2269674,
                  1740752, 1968596, 631872, 1734468,
                  1799540, 1910984, 398440, 2024770,
                  1940794, 1442408, 514826, 1403756,
                  2277398, 2318284, 774498, 2269674,
                  1740752, 1968596, 631872, 1734468,
                  1799540, 1910984, 398440, 2024770,
                  1299256, 1166912, 1555012, 1689928,
                  978400,  2119248, 2650754, 1684488,
                  884052,  1173252, 1993110, 1251508,
                  1391660, 1674368, 1559846, 2850934,
                  1299256, 1166912, 1555012, 1689928,
                  978400,  2119248, 2650754, 1684488,
                  884052,  1173252, 1993110, 1251508,
                  1391660, 1674368, 1559846, 2850934,
                  1299256, 1166912, 1555012, 1689928,
                  978400,  2119248, 2650754, 1684488,
                  884052,  1173252, 1993110, 1251508,
                  1391660, 1674368, 1559846, 2850934])
        else:
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
    #S = S[0:10,0:96]
    return S, index_signature


def read_counts(matrix_file):
    M = pd.read_csv(matrix_file, delimiter='\t')
    index_matrix = M.index.values.tolist()
    M = M.to_numpy().astype(float)
    #return pd.read_csv(matrix_file, delimiter='\t').to_numpy().astype(float)
    #print(index_matrix)
    return M, index_matrix



def get_num_cpus():
    return multiprocessing.cpu_count()


def denovo(matrix_file: str, n_signatures: int, lambd: float,
           opportunity_file: str = None, cosmic_file: str = None,
           max_em_iterations: int = 10000,
           max_gd_iterations: int = 50, numeric_chromosomes: bool = False, genotyped: bool = True, file_extension=None,
           ref_genome=None, output_folder: str = 'output/'):
    """
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
    numeric_chromosomes: bool
        True if chromosome names in vcf are '1', '2', '3'. False if 'chr1', 'chr2', 'chr3'
    genotyped: bool
        True if the VCF file has genotype information for many samples
    """
    start_time = time.time()
    num_cpus = get_num_cpus()
    print(f"Number of CPUs: {num_cpus}")
    if file_extension == '.vcf':
        assert ref_genome is not None, 'Please provide a reference genome along with the vcf file'
        count_mutation(matrix_file, ref_genome, f'{output_folder}/matrix.csv', numeric_chromosomes, genotyped)
        matrix_file = f'{output_folder}/matrix.csv'
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
        cos_similarity.to_csv(f"{output_folder}/cosine_similarity_denovo.txt", sep="\t")
        # print(cos_similarity)
    S = np.transpose(S)
    alphabet = list(string.ascii_uppercase)
    Sig = ['Denovo ' + alphabet[k] for k in range(S.shape[1])]
    # print(Sig)
    mutation_labels = ['A[C>A]A', 'A[C>A]C', 'A[C>A]G', 'A[C>A]T', 'C[C>A]A', 'C[C>A]C', 'C[C>A]G', 'C[C>A]T',
                       'G[C>A]A', 'G[C>A]C', 'G[C>A]G', 'G[C>A]T', 'T[C>A]A', 'T[C>A]C', 'T[C>A]G', 'T[C>A]T',
                       'A[C>G]A', 'A[C>G]C', 'A[C>G]G', 'A[C>G]T', 'C[C>G]A', 'C[C>G]C', 'C[C>G]G', 'C[C>G]T',
                       'G[C>G]A', 'G[C>G]C', 'G[C>G]G', 'G[C>G]T', 'T[C>G]A', 'T[C>G]C', 'T[C>G]G', 'T[C>G]T',
                       'A[C>T]A', 'A[C>T]C', 'A[C>T]G', 'A[C>T]T', 'C[C>T]A', 'C[C>T]C', 'C[C>T]G', 'C[C>T]T',
                       'G[C>T]A', 'G[C>T]C', 'G[C>T]G', 'G[C>T]T', 'T[C>T]A', 'T[C>T]C', 'T[C>T]G', 'T[C>T]T',
                       'A[T>A]A', 'A[T>A]C', 'A[T>A]G', 'A[T>A]T', 'C[T>A]A', 'C[T>A]C', 'C[T>A]G', 'C[T>A]T',
                       'G[T>A]A', 'G[T>A]C', 'G[T>A]G', 'G[T>A]T', 'T[T>A]A', 'T[T>A]C', 'T[T>A]G', 'T[T>A]T',
                       'A[T>C]A', 'A[T>C]C', 'A[T>C]G', 'A[T>C]T', 'C[T>C]A', 'C[T>C]C', 'C[T>C]G', 'C[T>C]T',
                       'G[T>C]A', 'G[T>C]C', 'G[T>C]G', 'G[T>C]T', 'T[T>C]A', 'T[T>C]C', 'T[T>C]G', 'T[T>C]T',
                       'A[T>G]A', 'A[T>G]C', 'A[T>G]G', 'A[T>G]T', 'C[T>G]A', 'C[T>G]C', 'C[T>G]G', 'C[T>G]T',
                       'G[T>G]A', 'G[T>G]C', 'G[T>G]G', 'G[T>G]T', 'T[T>G]A', 'T[T>G]C', 'T[T>G]G', 'T[T>G]T']
    label = list(mutation_labels)
    S = pd.DataFrame(S, columns=Sig, index=label)
    S.to_csv(f"{output_folder}/denovo_signature.txt", sep='\t', index=label)
    # print(S)
    deno_figure = plot_profile(S)
    deno_figure.savefig(f"{output_folder}/denovo_figure.png", dpi=600)
    E = pd.DataFrame(E, columns=Sig, index=None)
    # np.savetxt(output_file_exposure, np.array(E))
    E.to_csv(f"{output_folder}/denovo_exposures.txt", sep='\t', index=None)
    # np.savetxt(output_file_signature, np.array(S))
    print("--- %s seconds ---" % (time.time() - start_time))


def main():
    app = typer.Typer()
    app.command()(refit)
    app.command()(denovo)
    app.command()(count_mutation)
    app()


if __name__ == "__main__":
    # num_cpus = 4  # Define the number of CPUs to use
    # processes = []
    # for _ in range(num_cpus):
    #     process = multiprocessing.Process(target=denovo)
    #     processes.append(process)
    #     process.start()
    # for process in processes:
    #     process.join()
    main()
