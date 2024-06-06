"""Console script for mutational_starsign."""
from pathlib import Path
from typing import Optional
import logging
# from typing import Annotated
from typing_extensions import Annotated
import numpy as np
import pandas as pd
from enum import Enum
import time
import os
from matplotlib import pyplot as plt
import seaborn as sns
import string
import multiprocessing
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
np.random.seed(10000)

# todo


# class DataType(str, Enum):
#     exome = 'exome'
#     genome = 'genome'


import typer
from .refit import refit as _refit
from .denovo import denovo as _denovo, cos_sim_matrix
from .count_mutation_cli import count_mutation
from .bootstrapping import bootstrap as _bootstrap


def bootstrap(matrix_file: str, signature_file: str, output_file_exposure_avg: str, output_file_exposure_std: str,
              opportunity_file: str = None):
    M = read_counts(matrix_file)
    S, index_signature = read_signature(signature_file)
    O = read_opportunity(M, opportunity_file)
    ####    lambd = get_lambda(data_type)
    estimated_exposure, exposure_std = _bootstrap(M, S, O, lambd=lambd)
    np.savetxt(output_file_exposure_avg, estimated_exposure, delimiter='\t')
    np.savetxt(output_file_exposure_std, exposure_std, delimiter='\t')


def plot(file, output_folder='output/'):
    plt.style.use('default')
    file.plot(kind="bar", figsize=(10, 4))
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.ylabel("Mutation fraction")
    plt.xlabel("Signatures")
    plt.tight_layout()
    return plt.savefig(f"{output_folder}/plot_cohort.png", dpi=600)


def single_plot_old(file):  # 2 feb
    file = file.transpose()
    file.columns = ['Signatures', 'std_dev']
    # filtered_data = file[file['Signatures'] != 0]
    filtered_data = file[file['Signatures'] >= 0.05]
    color_palette = sns.color_palette("husl", n_colors=len(filtered_data))

    # Set up the figure and axis
    plt.style.use('default')
    fig, ax = plt.subplots()
    fig, ax = plt.subplots(figsize=(20, 8))

    # Plot the bar plot with error bars for non-zero values
    ax.bar(filtered_data.index, filtered_data['Signatures'], yerr=filtered_data['std_dev'],
           color=color_palette, alpha=0.5, align='center', capsize=3)
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

def single_plot_old2(file):
    data = file.transpose()
    filtered_data = data[data['Signature'] >= 0.06]
    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(20, 8))
    for i, row in filtered_data.iterrows():
        lower_error = row['E_25']
        upper_error = row['E_95']
        # lower_error = abs(row['Signature'] - row['E_25'])
        # upper_error = abs(row['E_95'] - row['Signature'])
        ax.bar(i, row['Signature'], yerr=[[lower_error], [upper_error]], capsize=5)
    ax.set_xticks(range(len(filtered_data)))
    ax.set_xticklabels(filtered_data.index)
    plt.xticks(rotation=45)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.xticks(rotation=45)
    ax.set_xlabel('Signature exposures')
    ax.set_ylabel('Mutation fraction')
    ax.set_title('Single Sample Mutational Signatures')
    return plt


def single_plot(data, target_median=0.06, show_percentiles=True):
    """
    This function creates a single violin plot for selected columns in a DataFrame,
    with column names on the x-axis and highlighting percentiles (optional).

    Args:
        data: A pandas DataFrame containing the data.
        target_median: The target median value for filtering columns (default: 0.06).
        show_percentiles: Boolean flag to display percentiles (default: True).

    Returns:
        A matplotlib figure object containing the violin plot.
    """

    # Calculate median for each column
    median_values = data.mean()

    # Filter columns based on median condition
    filtered_columns = median_values[median_values >= target_median].index.tolist()

    if not filtered_columns:
        print(f"No columns found with median value equal to {target_median}")
        return None

    # Subset the DataFrame to include only selected columns
    filtered_data = data[filtered_columns]

    # Reshape the DataFrame using melt to have all values in one column
    melted_data = filtered_data.melt()

    plt.figure(figsize=(12, 8))

    # Create violin plot
    ax = sns.violinplot(data=melted_data, x='variable', y='value', inner=None, cut=0)

    # Add swarmplot to show data points
    sns.swarmplot(data=melted_data, x='variable', y='value', color='k', size=3, ax=ax)

    # Add median as vertical bar inside violin
    for column_name in filtered_columns:
        median_val = filtered_data[column_name].median()
        ax.axvline(x=filtered_columns.index(column_name), ymin=0, ymax=1, color='white', linestyle='-', linewidth=2)
        #ax.axvline(x=filtered_columns.index(column_name), ymin=0.25, ymax=0.75, color='black', linestyle='-', linewidth=1)

    # Calculate and plot percentiles if needed
    if show_percentiles:
        percentiles = melted_data.groupby('variable')['value'].quantile([0.025, 0.975]).unstack()
        plt.scatter(range(len(filtered_columns)), percentiles.iloc[:, 0], color='red', marker='o', label='2.5th Percentile')
        plt.scatter(range(len(filtered_columns)), percentiles.iloc[:, 1], color='green', marker='o', label='97.5th Percentile')

    plt.xlabel('Signature exposures')
    plt.ylabel('Mutation fractions')
 #   plt.title('Violin Plot of Selected Columns')
    plt.xticks(ticks=range(len(filtered_columns)), labels=filtered_columns, rotation= 0, ha='right')
    if show_percentiles:
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    return plt

# Plot the violin plot
def single_plot_old(file):
    data_transposed = file.T
    data_transposed = data_transposed[data_transposed['Signature'] >= 0.06]
    data_transposed = data_transposed.T
    plt.figure(figsize=(10, 6))
    sns.violinplot(data=data_transposed, palette='Set2', cut=0)
    plt.title('Violin Plots of Signature with E_2.5 and E_97.5')
    plt.xlabel('Signature exposures')
    plt.ylabel('Mutation fraction')
    #plt.set_title('Single Sample Mutational Signatures')
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

    header = data.index
    data_list = data.values.tolist()
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

    fig = plt.figure(figsize=(15, 8))

    # Calculate number of rows correctly
    n_rows = len(data_list)
    fig_rows = int(np.ceil(n_rows / 2))  # Use ceiling to handle odd numbers

    for id_x in range(n_rows):
        plt.subplot(fig_rows, 2, id_x + 1)  # Use fig_rows directly
        plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.4)
        plt.bar(mutation_categories, data_list[id_x], color=mutation_colors)
        plt.xlabel('Mutation Categories')
        plt.ylabel('Mutations')
        plt.title(f'{header[id_x]}')
        plt.xticks(mutation_categories[::16], mutation_labels[::16], ha='center')
        plt.xlim(1, 97)
    #    plt.grid(b=True, color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)
        plt.grid(visible=True, color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)

        # Place legend outside of individual subplots
        if id_x == 0:
            legend_handles = [plt.Rectangle((0, 0), 1, 1, color=color) for color in color_groups.values()]
            legend_labels = color_groups.keys()
            plt.legend(legend_handles, legend_labels, bbox_to_anchor=(1.05, 1), loc='upper left')

    return plt






def refit(matrix_file: Annotated[str, typer.Argument(help='Tab separated matrix file')],
          signature_file: Annotated[str, typer.Argument(help='Tab separated matrix file')],
          ref_genome: str = None,
          n_bootstraps: int = 200,
          opportunity_file: Annotated[str, typer.Option(help="The distribution of triplets in a reference 'human-genome' or 'human-exome' or normal tissue")] = None,
          numeric_chromosomes: Annotated[bool, typer.Option(help="True if chromosome names in vcf are '1', '2', '3'. False if 'chr1', 'chr2', 'chr3'")] = False,
          genotyped: Annotated[bool, typer.Option(help="True if the VCF file has genotype information for many samples")] = True,
          output_folder: str = 'output/',
          signature_names: Annotated[Optional[str], typer.Option(help='Comma separated list of signature names')] = None,
          n_iterations: int=1000):

    '''
    Mutational Signatures Refit Parameters \n
    n_bootstraps (int): Number of bootstraps to consider for a single sample\n
    matrix_file (str): Path to a file containing Mutational Catalogue Matrix \n
    signature_file (str): Path to a file containing Known Mutational Signature reference signature\n
    output_file_exposure (str): Path to save the refitted exposure matrix \n
    opportunity_file (str): Path to a file defining Opportunity matrix\n
    numeric_chromosomes (bool): True if chromosome names in vcf are '1', '2', '3'. False if 'chr1', 'chr2', 'chr3' \n
    genotyped,bool (True) if the VCF file has genotype information for many samples \n
    ref_genome (str): Path to a file containing the reference genome\n
    n-iteration (int): Number of running iteration

    '''
    sig_name = Path(signature_file).stem
    matrix_name = Path(matrix_file).stem
    run_name = f'{matrix_name}_{sig_name}_{n_iterations}'
    logger.info(f'Refitting mutational signatures with {n_iterations} iterations')
    logger.info(f'Run name: {run_name}')
    start_time = time.time()
    # reference genome
    #ref_genome = '/Users/bope/Documents/MutSig/scientafellow/packages/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
    file_name, file_extension = os.path.splitext(matrix_file)
    if file_extension == '.vcf':
        assert ref_genome is not None, 'Please provide a reference genome along with the vcf file'
        count_mutation(matrix_file, ref_genome, f'{output_folder}/matrix.csv', numeric_chromosomes, genotyped)
        matrix_file = f'{output_folder}/matrix.csv'
    M = read_counts(matrix_file)
    index_matrix = M.index.values.tolist()
    S = read_signature(signature_file)
    # print(S)
    index_matrix = M.index.values.tolist()
    normalized_df = M.div(M.sum(axis=1), axis=0)
    col_means = normalized_df.mean(axis=0)
    count_less_than_01 = (col_means < 0.004).sum()
    threshold = 0.25 * len(col_means)
    if count_less_than_01 <= threshold:
        zero_contexts = np.array(M.columns)[(M.values < 0.01).any(axis=0)]
        # print("zero_contexts",zero_contexts)
        corr_sigs_mask = (S.loc[:, zero_contexts] >= 0.1).any(axis=1)
        signatures = S.loc[~S.index.isin(corr_sigs_mask[corr_sigs_mask].index)]
        S = signatures
    if signature_names is not None:
        S = filter_signatures(S, signature_names.split(','))

    index_signature = S.index.values.tolist()
    desired_order = M.columns

# Check for missing columns
    missing_cols = set(desired_order) - set(S.columns)  # Set difference for missing columns
    if len(missing_cols) > 0:
        raise ValueError("Error: Following columns are not present in the DataFrame:", missing_cols)
# Reorder columns
    S = S[desired_order]
    S = S.to_numpy().astype(float)
    M = M.to_numpy().astype(float)
    O = read_opportunity(M, opportunity_file)
    lambd = 0.7

    if (M.ndim == 2 and M.shape[0] == 1) or M.ndim == 1:
        boostrap_M = _bootstrap(M,n_bootstraps)

        expo_run = _refit(boostrap_M, S, O, lambd=lambd, n_iterations=n_iterations)
        expo_run = pd.DataFrame(data=expo_run, columns=index_signature)
        plot = single_plot(expo_run)
        expo_run_median = expo_run.median()
        expo_run_median.to_csv(f"{output_folder}/StarSign_exposure_median_{run_name}.txt", index=True, sep='\t')
        expo_run.to_csv(f"{output_folder}/StarSign_exposure_Exposure_test_{run_name}.txt", index=True, header=True, sep='\t')
        if plot is not None:
            plot.savefig(f"{output_folder}/StarSign_exposure_Exposure_test_{run_name}.png", dpi=600)
            plot.close()
        else:
            print("No plot was generated due to filtering criteria.")
    else:

        E = _refit(M, S, O, lambd=lambd, n_iterations=n_iterations)
        # print(O)
        sum_expo = E.sum(axis=0, keepdims=True) / len(E)
        sum_expo_t = np.transpose(sum_expo)
        E = pd.DataFrame(data=E, columns=index_signature, index=index_matrix)
        E.to_csv(f'{output_folder}/{run_name}.txt', index=index_matrix, header=True, sep='\t')
        sum_expo = pd.DataFrame(data=sum_expo, columns=index_signature, index=['Signatures'])
        sum_expo.to_csv(f'{output_folder}/average_{run_name}.txt',columns=index_signature, index= False, sep='\t')
        sum_expo = np.transpose(sum_expo)
        sort_E = sum_expo.sort_values(by=['Signatures'], ascending=False)
        sort_E = sort_E.iloc[:5, 0:]
        index_sort_E = sort_E.index.values.tolist()
        sort_E = pd.DataFrame(data=sort_E, index=index_sort_E)
        sort_E = sort_E.T
        plot_top_five = cohort_violin(sort_E)
        plot_top_five.savefig(f"{output_folder}/starsign_top5_signatures_{run_name}.png", dpi=600)
        plot_variance = cohort_violin(E)
        plot_variance.savefig(f"{output_folder}/starsign_cohort_{run_name}.png", dpi=600)
    print("--- %s seconds ---" % (time.time() - start_time))



def read_opportunity(M, opportunity_file):
    n_samples = len(M)
    n_mutations = M.shape[1]
    if opportunity_file is not None:
        if opportunity_file == "human-genome":
            O = np.array([1.14e+08, 6.60e+07, 1.43e+07, 9.12e+07,  # C>A @ AC[ACGT]
                          1.05e+08, 7.46e+07, 1.57e+07, 1.01e+08,  # C>A @ CC[ACGT]
                          8.17e+07, 6.76e+07, 1.35e+07, 7.93e+07,  # C>A @ GC[ACGT]
                          1.11e+08, 8.75e+07, 1.25e+07, 1.25e+08,  # C>A @ TC[ACGT]
                          1.14e+08, 6.60e+07, 1.43e+07, 9.12e+07,  # C>G @ AC[ACGT]
                          1.05e+08, 7.46e+07, 1.57e+07, 1.01e+08,  # C>G @ CC[ACGT]
                          8.17e+07, 6.76e+07, 1.35e+07, 7.93e+07,  # C>G @ GC[ACGT]
                          1.11e+08, 8.75e+07, 1.25e+07, 1.25e+08,  # C>G @ TC[ACGT]
                          1.14e+08, 6.60e+07, 1.43e+07, 9.12e+07,  # C>T @ AC[ACGT]
                          1.05e+08, 7.46e+07, 1.57e+07, 1.01e+08,  # C>T @ CC[ACGT]
                          8.17e+07, 6.76e+07, 1.35e+07, 7.93e+07,  # C>T @ GC[ACGT]
                          1.11e+08, 8.75e+07, 1.25e+07, 1.25e+08,  # C>T @ TC[ACGT]
                          1.17e+08, 7.57e+07, 1.04e+08, 1.41e+08,  # T>A @ AC[ACGT]
                          7.31e+07, 9.55e+07, 1.15e+08, 1.13e+08,  # T>A @ CC[ACGT]
                          6.43e+07, 5.36e+07, 8.52e+07, 8.27e+07,  # T>A @ GC[ACGT]
                          1.18e+08, 1.12e+08, 1.07e+08, 2.18e+08,  # T>A @ TC[ACGT]
                          1.17e+08, 7.57e+07, 1.04e+08, 1.41e+08,  # T>C @ AC[ACGT]
                          7.31e+07, 9.55e+07, 1.15e+08, 1.13e+08,  # T>C @ CC[ACGT]
                          6.43e+07, 5.36e+07, 8.52e+07, 8.27e+07,  # T>C @ GC[ACGT]
                          1.18e+08, 1.12e+08, 1.07e+08, 2.18e+08,  # T>C @ TC[ACGT]
                          1.17e+08, 7.57e+07, 1.04e+08, 1.41e+08,  # T>G @ AC[ACGT]
                          7.31e+07, 9.55e+07, 1.15e+08, 1.13e+08,  # T>G @ AC[ACGT]
                          6.43e+07, 5.36e+07, 8.52e+07, 8.27e+07,  # T>G @ AG[ACGT]
                          1.18e+08, 1.12e+08, 1.07e+08, 2.18e+08])  # T>G @ AT[ACGT]])
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
                          978400, 2119248, 2650754, 1684488,
                          884052, 1173252, 1993110, 1251508,
                          1391660, 1674368, 1559846, 2850934,
                          1299256, 1166912, 1555012, 1689928,
                          978400, 2119248, 2650754, 1684488,
                          884052, 1173252, 1993110, 1251508,
                          1391660, 1674368, 1559846, 2850934,
                          1299256, 1166912, 1555012, 1689928,
                          978400, 2119248, 2650754, 1684488,
                          884052, 1173252, 1993110, 1251508,
                          1391660, 1674368, 1559846, 2850934])
        else:
            O = pd.read_csv(opportunity_file, sep='\t', header=None).to_numpy().astype(float)
            normalized_vector1 = O / np.linalg.norm(O)
            min_value_vector2 = np.min(M)
            max_value_vector2 = np.max(M)
            O = normalized_vector1 * (max_value_vector2 - min_value_vector2) + min_value_vector2
        O = np.broadcast_to(O, M.shape)
    else:
        O = np.ones((n_samples, n_mutations), dtype=float)
    O = O / np.amin(O).sum(axis=-1, keepdims=True)

    assert O.shape == (n_samples, n_mutations), f'{O.shape} != {(n_samples, n_mutations)}'
    return O


def read_signature(signature_file):
    S = pd.read_csv(signature_file, delimiter='\t')
    return S


def read_counts(matrix_file):
    M = pd.read_csv(matrix_file, delimiter='\t')
    return M


def filter_signatures(S, signature_names):
    if len(signature_names) < 3:
        raise ValueError("You must select at least 3 signature names.")
    S = S.loc[signature_names]
    return S


def get_tri_context_fraction(mut_counts):
    total_mutations = mut_counts.sum().sum()
    trinucleotide_fractions = mut_counts.div(total_mutations) if total_mutations != 0 else mut_counts
    return trinucleotide_fractions


def get_num_cpus():
    return multiprocessing.cpu_count()


def denovo(matrix_file: Annotated[str, typer.Argument(help='Tab separated matrix file')],
           n_signatures: Annotated[int, typer.Argument(help='Tab separated signature file')],
           lambd: Annotated[float, typer.Option(help='Regularization parameter')] = 0.7,
       #    opportunity_file: str = None,
           opportunity_file: Annotated[str, typer.Option(help="The distribution of triplets in a reference 'human-genome' or 'human-exome' or normal tissue")] = None,
         #  opportunity_file: Annotated[str,typer.Argument(help='The distribution of triplets in a reference genome/exome or normal tissue')] = None,
           cosmic_file: Annotated[str, typer.Option(help='Tab separated cosmic file')] = None,
           numeric_chromosomes: Annotated[bool, typer.Option(help="True if chromosome names in vcf are '1', '2', '3'. False if 'chr1', 'chr2', 'chr3'")] = False,
           genotyped: Annotated[bool, typer.Option(help="True if the VCF file has genotype information for many samples")] = True,
           max_em_iterations: int = 100,
           max_gd_iterations: int = 50,
           file_extension=None,
           ref_genome=None, output_folder: str = 'output/'):

    """Performs denovo  Mutational Signatures analysis.

    Args:   matrix_file (str): Path to the tab-separated matrix file \n
            n_signatures (int): Number of signatures to identify \n
            lambd (float, optional): Regularization parameter. Defaults to 0.7 \n
            help_lambda (bool, optional): Flag to display lambda help message. Defaults to False \n
            numeric_chromosomes (bool): True if chromosome names in vcf are '1', '2', '3'. False if 'chr1', 'chr2', 'chr3' \n
            genotyped (bool) : True if the VCF file has genotype information for many samples
    """
    matrix_name = Path(matrix_file).stem
    run_name = f'{matrix_name}'
    logger.info(f'Run name: {run_name}')
    start_time = time.time()
    if file_extension == '.vcf':
        assert ref_genome is not None, 'Please provide a reference genome along with the vcf file'
        count_mutation(matrix_file, ref_genome, f'{output_folder}/matrix.csv', numeric_chromosomes, genotyped)
        matrix_file = f'{output_folder}/matrix.csv'
    M  = read_counts(matrix_file)
    index_matrix = M.index.values.tolist()
    n_samples = len(M)
    n_signatures = n_signatures
    lambd = lambd
    print('lambda',lambd)
    n_mutations = M.shape[1]
    O = read_opportunity(M, opportunity_file)
    desired_order = M.columns
    #print(desired_order)
    M = M.to_numpy().astype(float)
    E, S = _denovo(M, n_signatures, lambd, O, em_steps=max_em_iterations, gd_steps=max_gd_iterations)
    if cosmic_file is not None:
        cosmic = pd.read_csv(cosmic_file, delimiter='\t')

# Check for missing columns
        missing_cols = set(desired_order) - set(cosmic.columns)  # Set difference for missing columns
        if len(missing_cols) > 0:
            raise ValueError("Error: Following columns are not present in the DataFrame:", missing_cols)
# Reorder columns
        cosmic = cosmic[desired_order]
        cos_similarity = cos_sim_matrix(S, cosmic)[0]
        cos_similarity.to_csv(f"{output_folder}/StarSign_{run_name}_cosine_similarity.txt", sep="\t")
        # print(cos_similarity)
    S = np.transpose(S)
    alphabet = list(string.ascii_uppercase)
    Sig = ['De novo ' + alphabet[k] for k in range(S.shape[1])]
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
    S.to_csv(f"{output_folder}/StarSign_{run_name}_Denovo_signature.txt", sep='\t', index=label)
    deno_figure = plot_profile(S.T)
    deno_figure.savefig(f"{output_folder}/StarSign_{run_name}_profile.png", dpi=600)
    E = pd.DataFrame(E, columns=Sig, index=None)
    E.to_csv(f"{output_folder}/StarSign_{run_name}_Denovo_exposures.txt", sep='\t', index=None)
    print("--- %s seconds ---" % (time.time() - start_time))


def main():
    app = typer.Typer()
    app.command()(refit)
    app.command()(denovo)
    app.command()(count_mutation)
    app()


if __name__ == "__main__":
    main()
