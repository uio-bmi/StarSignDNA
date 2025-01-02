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
np.random.seed(1000)

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

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


def plot(file):
    plt.style.use('default')
    file.plot(kind="bar", figsize=(10, 4))
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.ylabel("Mutation fraction")
    plt.xlabel("Signatures")
    plt.tight_layout()
    return plt



def single_plot(data):
    df = pd.DataFrame(data)

    # Calculate medians
    medians = df.median()

    # Select columns where the median is above the threshold
    threshold = 0.06
    selected_columns = medians[medians > threshold].index

    # Filter the DataFrame
    filtered_df = df[selected_columns]
# Calculate percentiles for the filtered DataFrame
    percentiles = filtered_df.quantile([0.025, 0.975])

    # Create the plot
    plt.figure(figsize=(12, 8))

    # Plotting the data points
    for column in filtered_df.columns:
        plt.violinplot(filtered_df[column], positions=[list(filtered_df.columns).index(column)], showmeans=True)
        plt.scatter([list(filtered_df.columns).index(column)], [percentiles.at[0.025, column]], color='red', zorder=5, label='2.5th Percentile' if column == filtered_df.columns[0] else "")
        plt.scatter([list(filtered_df.columns).index(column)], [percentiles.at[0.975, column]], color='green', zorder=5, label='97.5th Percentile' if column == filtered_df.columns[0] else "")

    plt.xticks(ticks=range(len(filtered_df.columns)), labels=filtered_df.columns)
    plt.xlabel('Signature Exposures')
    plt.ylabel('Mutation Fractions')
    plt.legend()
    return plt
# Plot the violin plot



# Version with more customization options
def single_plot_custom(data, figsize=(12, 8),
                      violin_color='#80B1D3',
                      violin_alpha=0.3,
                      median_color='#80B1D3',
                      error_bar_color='black',
                      error_bar_width=2):
    df = pd.DataFrame(data)
    medians = df.median()

    # Filter based on threshold
    threshold = 0.1
    selected_columns = medians[medians > threshold].index
    filtered_df = df[selected_columns]

    # Calculate statistics
    percentiles = filtered_df.quantile([0.025, 0.975])
    medians = filtered_df.median()

    # Create plot
    fig, ax = plt.subplots(figsize=figsize)

    # Plot violins and error bars
    for i, column in enumerate(filtered_df.columns):
        # Violin plot
        parts = ax.violinplot(filtered_df[column],
                            positions=[i],
                            showmeans=False,
                            showextrema=False)

        # Customize violin
        for pc in parts['bodies']:
            pc.set_facecolor(violin_color)
            pc.set_alpha(violin_alpha)

        # Error bars
        ax.vlines(x=i,
                 ymin=percentiles.at[0.025, column],
                 ymax=percentiles.at[0.975, column],
                 color=error_bar_color,
                 linewidth=error_bar_width,
                 zorder=3)

        # Median points
        ax.plot(i, medians[column], 'o',
               color=median_color,
               markersize=8,
               zorder=4)

    # Customize plot
    ax.set_xticks(range(len(filtered_df.columns)))
    ax.set_xticklabels(filtered_df.columns, rotation=45, ha='right')
    ax.set_xlabel('Signature Exposures')
    ax.set_ylabel('Mutation Fractions')

    # Add grid
    ax.grid(True, linestyle='--', alpha=0.3)

    # Add legend
    ax.vlines([], [], [], color=error_bar_color, linewidth=error_bar_width, label='95% CI')
    ax.plot([], [], 'o', color=median_color, label='Median')
    ax.legend()

    plt.tight_layout()

    return plt

# Optional: Add more customization
def single_plot_custom(data, figsize=(12, 8), median_color='#80B1D3',
                      error_bar_color='black', error_bar_width=2):
    df = pd.DataFrame(data)
    medians = df.median()

    # Filter based on threshold
    threshold = 0.1
    selected_columns = medians[medians > threshold].index
    filtered_df = df[selected_columns]

    # Calculate statistics
    percentiles = filtered_df.quantile([0.025, 0.975])
    medians = filtered_df.median()

    # Create plot
    fig, ax = plt.subplots(figsize=figsize)

    # Plot error bars and medians
    for i, column in enumerate(filtered_df.columns):
        # Error bars
        ax.vlines(x=i,
                 ymin=percentiles.at[0.025, column],
                 ymax=percentiles.at[0.975, column],
                 color=error_bar_color,
                 linewidth=error_bar_width)

        # Median points
        ax.plot(i, medians[column], 'o',
               color=median_color,
               markersize=8,
               zorder=5)

    # Customize plot
    ax.set_xticks(range(len(filtered_df.columns)))
    ax.set_xticklabels(filtered_df.columns, rotation=45, ha='right')
    ax.set_xlabel('Signature Exposures')
    ax.set_ylabel('Mutation Fractions')

    # Add grid
    ax.grid(True, linestyle='--', alpha=0.3)

    # Add confidence interval explanation
    ax.text(0.02, 0.98, '95% Confidence Interval',
            transform=ax.transAxes,
            verticalalignment='top',
            fontsize=10)

    plt.tight_layout()

    return plt

def cohort_plot(file):
    # Get the colorblind-friendly palette from seaborn
    num_colors = len(file.index)
    color_palette = sns.color_palette("colorblind", num_colors)
    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(10, 4))
    # Transpose the dataframe so that the signatures are columns
    file.plot(kind="bar", color=color_palette, ax=ax)
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


def download_reference_genome(ref_genome=None, genome_path=None, dest_dir='genomes'):
    import os
    import requests
    import gzip
    import shutil
    import typer
    from typing import Optional, Annotated
    """
    Download and unzip a reference genome if not already available locally,
    or use a provided reference genome path.

    Parameters:
    ref_genome (str): The name of the reference genome to download (e.g., 'GRCh37', 'mm10').
    genome_path (str): The local path to the reference genome file.
    dest_dir (str): The directory to save the downloaded genome file. Default is 'genomes'.

    Returns:
    str: The path to the reference genome file.
    """
    if genome_path:
        if os.path.isfile(genome_path):
            print(f"Using provided reference genome at {genome_path}")
            return genome_path
        else:
            raise FileNotFoundError(f"The provided genome path {genome_path} does not exist.")

    if ref_genome is None:
        raise ValueError("Either ref_genome or genome_path must be provided.")

    genome_urls = {
        'GRCh37': 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz',
        'GRCh38': 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz',
        'mm10': 'http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz',
        'mm9': 'http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.fa.gz',
        'rn6': 'http://hgdownload.cse.ucsc.edu/goldenPath/rn6/bigZips/rn6.fa.gz',
        'mm6': 'http://hgdownload.cse.ucsc.edu/goldenPath/mm6/bigZips/mm6.fa.gz'
    }

    if ref_genome not in genome_urls:
        raise ValueError(f"Reference genome {ref_genome} is not supported.")

    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)

    genome_file_gz = os.path.join(dest_dir, f"{ref_genome}.fa.gz")
    genome_file = os.path.join(dest_dir, f"{ref_genome}.fa")

    if not os.path.isfile(genome_file):
        if not os.path.isfile(genome_file_gz):
            print(f"Downloading reference genome {ref_genome}...")
            response = requests.get(genome_urls[ref_genome], stream=True)
            if response.status_code == 200:
                with open(genome_file_gz, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                print(f"Downloaded {ref_genome} successfully.")
            else:
                raise Exception(f"Failed to download reference genome {ref_genome}.")
        else:
            print(f"Compressed reference genome {ref_genome} already exists locally.")

        print(f"Unzipping {genome_file_gz}...")
        with gzip.open(genome_file_gz, 'rb') as f_in:
            with open(genome_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print(f"Unzipped {genome_file} successfully.")
    else:
        print(f"Reference genome {ref_genome} already exists locally.")

    return genome_file


def refit(matrix_file: Annotated[str, typer.Argument(help='Tab separated matrix file')],
          signature_file: Annotated[str, typer.Argument(help='Tab separated matrix file')],
          ref_genome: Annotated[Optional[str], typer.Option(help='Reference genome to use. Options: GRCh37, GRCh38, mm10, mm9, rn6, mm6')] = None,
          genome_path: Annotated[Optional[str], typer.Option(help='Local path to the reference genome file')] = None,
          n_bootstraps: int = 200,
          opportunity_file: Annotated[str, typer.Option(help="The distribution of triplets in a reference 'human-genome' or 'human-exome' or normal tissue")] = None,
          numeric_chromosomes: Annotated[bool, typer.Option(help="True if chromosome names in vcf are '1', '2', '3'. False if 'chr1', 'chr2', 'chr3'")] = True,
          genotyped: Annotated[bool, typer.Option(help="True if the VCF file has genotype information for many samples")] = True,
          output_folder: str = 'output/',
          signature_names: Annotated[Optional[str], typer.Option(help='Comma separated list of signature names')] = None,
          n_iterations: int= 1000,
          data_type: DataType = DataType.exome):


    '''
    Mutational Signatures Refit Parameters \n
    n_bootstraps (int): Number of bootstraps to consider for a single sample\n
    matrix_file (str): Path to a file containing Mutational Catalogue Matrix \n
    signature_file (str): Path to a file containing Known Mutational Signature reference signature\n
    output_file_exposure (str): Path to save the refitted exposure matrix \n
    opportunity_file (str): Path to a file defining weight matrix\n
    numeric_chromosomes (bool): True if chromosome names in vcf are '1', '2', '3'. False if 'chr1', 'chr2', 'chr3' \n
    genotyped,bool (True) if the VCF file has genotype information for many samples \n
    ref_genome (str): Path to a file containing the reference genome\n
    n-iteration (int): Number of running iteration
    data_type: The default lambda values are set as follows: Exome: 0.7, Genome: 100. These values are optimized for the respective data types to ensure accurate analysis and results.

    '''
    sig_name = Path(signature_file).stem
    matrix_name = Path(matrix_file).stem
    run_name = f'{matrix_name}_{sig_name}_{n_iterations}'
    logger.info(f'Refitting mutational signatures with {n_iterations} iterations')
    logger.info(f'Run name: {run_name}')
    start_time = time.time()
    file_name, file_extension = os.path.splitext(matrix_file)
    provided_genome_path = None
    if file_extension == '.vcf':
        if not ref_genome and not genome_path:
            raise ValueError("Either ref_genome or genome_path must be provided.")
        genome_path = download_reference_genome(ref_genome=ref_genome, genome_path=genome_path)
        print(f"Reference genome path: {genome_path}")
        count_mutation(matrix_file, genome_path, f'{output_folder}/matrix.csv', numeric_chromosomes, genotyped)
        matrix_file = f'{output_folder}/matrix.csv'
    M = read_counts(matrix_file)
    index_matrix = M.index.values.tolist()
    S = read_signature(signature_file)
    # print(S)

    #comment
    index_matrix = M.index.values.tolist()
    normalized_df = M.div(M.sum(axis=1), axis=0)
    col_means = normalized_df.mean(axis=0)
    count_less_than_01 = (col_means < 0.004).sum()
    threshold = 0.25 * len(col_means)
    if count_less_than_01 <= threshold:
        zero_contexts = np.array(M.columns)[(M.values < 0.01).any(axis=0)]
        # print("zero_contexts",zero_contexts)
        corr_sigs_mask = (S.loc[:, zero_contexts] >= 0.06).any(axis=1)
        signatures = S.loc[~S.index.isin(corr_sigs_mask[corr_sigs_mask].index)]
        S = signatures
    if signature_names is not None:
        S = filter_signatures(S, signature_names.split(','))
#comment
    index_signature = S.index.values.tolist()
    desired_order = M.columns

# Check for missing columns
    missing_cols = set(desired_order) - set(S.columns)  # Set difference for missing columns
    if len(missing_cols) > 0:
        raise ValueError("Error: Following columns are not present in the DataFrame or make sure the file provided is Tab separated:", missing_cols)
# Reorder columns
    S = S[desired_order]
    S = S.to_numpy().astype(float)
    M = M.to_numpy().astype(float)
    O = read_opportunity(M, opportunity_file)
    lambd = get_lambda(data_type)
    if (M.ndim == 2 and M.shape[0] == 1) or M.ndim == 1:
        boostrap_M = _bootstrap(M,n_bootstraps)
        print("Lambda",lambd)

        expo_run = _refit(boostrap_M, S, O, lambd=lambd, n_iterations=n_iterations)
        expo_run = pd.DataFrame(data=expo_run, columns=index_signature)
        #print(expo_run)

        mean_columns_E = expo_run.mean(axis=0)

        for row in M:
            row_sum = row.sum()
            if row_sum > 600:
                threshold = 0.02
            else:
                threshold = 0.03
       # threshold = 0.02
        selected_headers = mean_columns_E[mean_columns_E >= threshold].index.tolist()
        S = read_signature(signature_file)
        S = S.loc[selected_headers]
        index_signature = S.index.values.tolist()
        expo_run = _refit(boostrap_M, S, O, lambd=lambd, n_iterations=n_iterations)
        expo_run = pd.DataFrame(data=expo_run, columns=index_signature)
        print(expo_run)
        plot = single_plot(expo_run)
        expo_run_median = expo_run.median()
        expo_run_median.to_csv(f"{output_folder}/StarSign_exposure_median_{run_name}.txt", index=True, sep='\t')
        expo_run.to_csv(f"{output_folder}/StarSign_exposure_Exposure_{run_name}.txt", index=True, header=True, sep='\t')
        if plot is not None:
            plot.savefig(f"{output_folder}/StarSign_exposure_Exposure_{run_name}.png", dpi=300)
            plot.close()
        else:
            print("No plot was generated due to filtering criteria.")
    else:
        print("original",S.shape)
        E = _refit(M, S, O, lambd=lambd, n_iterations=n_iterations)
        # print(O)

        E = pd.DataFrame(data=E, columns=index_signature, index=index_matrix)
        print(E)
        # 1. Compute the mean of each column in E
        mean_columns_E = E.mean(axis=0)

        for row in M:
            row_sum = row.sum()
            if row_sum > 600:
                threshold = 0.02
            else:
                threshold = 0.03
       # threshold = 0.02
        selected_headers = mean_columns_E[mean_columns_E >= threshold].index.tolist()
        S = read_signature(signature_file)
        S = S.loc[selected_headers]

        index_signature = S.index.values.tolist()
        E = _refit(M, S, O, lambd=lambd, n_iterations=n_iterations)
        sum_expo = E.sum(axis=0, keepdims=True) / len(E)
        E = pd.DataFrame(data=E, columns=index_signature, index=index_matrix)
        E.to_csv(f'{output_folder}/{run_name}_threshold.txt', index=index_matrix, header=True, sep='\t')
        sum_expo = pd.DataFrame(data=sum_expo, columns=index_signature, index=['Signatures'])
        sum_expo.to_csv(f'{output_folder}/average_{run_name}.txt',columns=index_signature, index= False, sep='\t')
        sum_expo = np.transpose(sum_expo)
        sort_E = sum_expo.sort_values(by=['Signatures'], ascending=False)
        sort_E = sort_E.iloc[:5, 0:]
        index_sort_E = sort_E.index.values.tolist()
        sort_E = pd.DataFrame(data=sort_E, index=index_sort_E)
        plot_top_five = cohort_plot(sort_E)
        plot_top_five.savefig(f"{output_folder}/starsign_top5_signatures_{run_name}.png", dpi=300)
        plot_variance = cohort_violin(E)
        plot_variance.savefig(f"{output_folder}/starsign_cohort_{run_name}.png", dpi=300)
    print("--- %s seconds ---" % (time.time() - start_time))


def get_lambda(data_type):
    if data_type == DataType.genome:
        lambd = 1000
    else:
        lambd = 0.7
    return lambd

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

            def zscore_normalize(df):
                return (df - df.mean()) / df.std()

            def context_frequency_normalization(mutation_matrix, context_frequencies):
  # Convert context_frequencies to an array to divide by each column in mutation_matrix
                context_freq_array = np.array([context_frequencies[key] for key in context_frequencies.keys()])

# Normalize each column in the mutation_matrix by the corresponding context frequency
                normalized_matrix = mutation_matrix / context_freq_array
            return O
            context_frequencies = {'ACA>C': 59474352,'CGA>T': 6541004,'TGT>A': 59845313,'GCT>C': 40767054
}
            O = context_frequency_normalization(O, context_frequencies)
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
    if len(signature_names) < 5:
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
           ref_genome: Annotated[Optional[str], typer.Option(help='Reference genome to use. Options: GRCh37, GRCh38, mm10, mm9, rn6, mm6')] = None,
           genome_path: Annotated[Optional[str], typer.Option(help='Local path to the reference genome file')] = None,
           opportunity_file: Annotated[str, typer.Option(help="The distribution of triplets in a reference 'human-genome' or 'human-exome' or normal tissue")] = None,
           cosmic_file: Annotated[str, typer.Option(help='Tab separated cosmic file')] = None,
           numeric_chromosomes: Annotated[bool, typer.Option(help="True if chromosome names in vcf are '1', '2', '3'. False if 'chr1', 'chr2', 'chr3'")] = False,
           genotyped: Annotated[bool, typer.Option(help="True if the VCF file has genotype information for many samples")] = True,
           max_em_iterations: int = 100,
           max_gd_iterations: int = 50,
           file_extension=None,
           output_folder: str = 'output/'):

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
        if not ref_genome and not genome_path:
            raise ValueError("Either ref_genome or genome_path must be provided.")
        genome_path = download_reference_genome(ref_genome=ref_genome, genome_path=genome_path)
        print(f"Reference genome path: {genome_path}")
        count_mutation(matrix_file, genome_path, f'{output_folder}/matrix.csv', numeric_chromosomes, genotyped)
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
