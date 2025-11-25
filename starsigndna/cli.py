"""Command line interface for mutational signature analysis."""
from pathlib import Path
from typing import Optional
import logging
from typing_extensions import Annotated
import numpy as np
import pandas as pd
from enum import Enum
import time
import os
import matplotlib.pyplot as plt
import seaborn as sns
import string
import multiprocessing
import requests
import gzip
import shutil
import typer

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
np.random.seed(1000)

from .refit import refit as _refit
from .denovo import denovo as _denovo, cos_sim_matrix
from .count_mutation_cli import count_mutation
from .bootstrapping import bootstrap as _bootstrap

class DataType(str, Enum):
    """Enum for data types used in analysis."""
    exome = 'exome'
    genome = 'genome'

def plot(file):
    """Create a basic bar plot of mutation fractions.
    
    Args:
        file: DataFrame containing mutation data
    Returns:
        matplotlib.pyplot figure
    """
    plt.style.use('default')
    file.plot(kind="bar", figsize=(10, 4))
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.ylabel("Mutation fraction")
    plt.xlabel("Signatures")
    plt.tight_layout()
    return plt

def single_plot(data):
    """Create violin plot for signature exposures.
    
    Args:
        data: DataFrame containing signature exposure data
    Returns:
        matplotlib.pyplot figure
    """
    df = pd.DataFrame(data)
    medians = df.median()
    
    # Filter based on threshold
    threshold = 0.06
    selected_columns = medians[medians > threshold].index
    filtered_df = df[selected_columns]
    
    # Calculate percentiles
    percentiles = filtered_df.quantile([0.025, 0.975])
    
    plt.figure(figsize=(12, 8))
    for column in filtered_df.columns:
        plt.violinplot(filtered_df[column], positions=[list(filtered_df.columns).index(column)], showmeans=True)
        plt.scatter([list(filtered_df.columns).index(column)], [percentiles.at[0.025, column]], color='red', zorder=5, 
                   label='2.5th Percentile' if column == filtered_df.columns[0] else "")
        plt.scatter([list(filtered_df.columns).index(column)], [percentiles.at[0.975, column]], color='green', zorder=5, 
                   label='97.5th Percentile' if column == filtered_df.columns[0] else "")

    plt.xticks(ticks=range(len(filtered_df.columns)), labels=filtered_df.columns)
    plt.xlabel('Signature Exposures')
    plt.ylabel('Mutation Fractions')
    plt.legend()
    return plt

def single_plot_custom(data, figsize=(12, 8), violin_color='#80B1D3', violin_alpha=0.3,
                      median_color='#80B1D3', error_bar_color='black', error_bar_width=2):
    """Create customized violin plot for signature exposures.
    
    Args:
        data: DataFrame containing signature exposure data
        figsize: Figure size tuple
        violin_color: Color for violin plots
        violin_alpha: Transparency for violin plots
        median_color: Color for median points
        error_bar_color: Color for error bars
        error_bar_width: Width of error bars
    Returns:
        matplotlib.pyplot figure
    """
    df = pd.DataFrame(data)
    medians = df.median()
    
    # Filter based on threshold
    threshold = 0.1
    selected_columns = medians[medians > threshold].index
    filtered_df = df[selected_columns]
    
    # Calculate statistics
    percentiles = filtered_df.quantile([0.025, 0.975])
    medians = filtered_df.median()
    
    fig, ax = plt.subplots(figsize=figsize)
    
    for i, column in enumerate(filtered_df.columns):
        # Error bars
        ax.vlines(x=i, ymin=percentiles.at[0.025, column],
                 ymax=percentiles.at[0.975, column],
                 color=error_bar_color, linewidth=error_bar_width)
        
        # Median points
        ax.plot(i, medians[column], 'o', color=median_color,
               markersize=8, zorder=5)
    
    ax.set_xticks(range(len(filtered_df.columns)))
    ax.set_xticklabels(filtered_df.columns, rotation=45, ha='right')
    ax.set_xlabel('Signature Exposures')
    ax.set_ylabel('Mutation Fractions')
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.text(0.02, 0.98, '95% Confidence Interval',
            transform=ax.transAxes, verticalalignment='top', fontsize=10)
    
    plt.tight_layout()
    return plt

def cohort_plot(file):
    """Create bar plot for cohort analysis using colorblind-friendly palette.
    
    Args:
        file: DataFrame containing cohort data
    Returns:
        matplotlib.pyplot figure or None if data is empty
    """
    # Check if DataFrame is empty
    if file.empty or len(file.index) == 0:
        logger.warning("Cannot create cohort plot: DataFrame is empty")
        return None
    
    num_colors = len(file.index)
    color_palette = sns.color_palette("colorblind", num_colors)
    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(10, 4))
    file.plot(kind="bar", color=color_palette, ax=ax)
    plt.xticks(rotation=45, fontsize=7)
    plt.yticks(fontsize=7)
    plt.ylabel("Mutation fraction")
    plt.xlabel("Signature exposures")
    plt.tight_layout()
    return plt

def cohort_violin(file):
    """Create violin plot for cohort analysis.
    
    Args:
        file: DataFrame containing cohort data
    Returns:
        matplotlib.pyplot figure or None if data is empty
    """
    # Check if DataFrame is empty
    if file.empty or len(file.columns) == 0:
        logger.warning("Cannot create cohort violin plot: DataFrame is empty")
        return None
    
    fig, ax = plt.subplots(figsize=(20, 8))
    sns.violinplot(data=file, ax=ax, scale="count", cut=0)
    plt.xticks(rotation=45, fontsize=8)
    plt.yticks(fontsize=8)
    plt.ylabel("Mutation fraction")
    plt.xlabel("Signature exposures")
    return plt

def plot_profile(data):
    """Create mutation profile plot.
    
    Args:
        data: DataFrame containing mutation profile data
    Returns:
        matplotlib.pyplot figure
    """
    plt.style.use('default')
    header = data.index
    data_list = data.values.tolist()
    mutation_categories = np.arange(1, 97)
    
    # Define mutation labels and colors
    mutation_labels = ['C>A', 'C>A', 'C>A', 'C>A', 'C>A', 'C>A', 'C>A', 'C>A',
                      'C>G', 'C>G', 'C>G', 'C>G', 'C>G', 'C>G', 'C>G', 'C>G',
                      'C>T', 'C>T', 'C>T', 'C>T', 'C>T', 'C>T', 'C>T', 'C>T',
                      'T>A', 'T>A', 'T>A', 'T>A', 'T>A', 'T>A', 'T>A', 'T>A',
                      'T>C', 'T>C', 'T>C', 'T>C', 'T>C', 'T>C', 'T>C', 'T>C',
                      'T>G', 'T>G', 'T>G', 'T>G', 'T>G', 'T>G', 'T>G', 'T>G'] * 2

    color_groups = {
        'C>A': 'red', 'C>G': 'blue', 'C>T': 'green',
        'T>A': 'orange', 'T>C': 'purple', 'T>G': 'brown'
    }
    mutation_colors = [color_groups[label] for label in mutation_labels]

    fig = plt.figure(figsize=(15, 8))
    n_rows = len(data_list)
    fig_rows = int(np.ceil(n_rows / 2))

    for id_x in range(n_rows):
        plt.subplot(fig_rows, 2, id_x + 1)
        plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.4)
        plt.bar(mutation_categories, data_list[id_x], color=mutation_colors)
        plt.xlabel('Mutation Categories')
        plt.ylabel('Mutations')
        plt.title(f'{header[id_x]}')
        plt.xticks(mutation_categories[::16], mutation_labels[::16], ha='center')
        plt.xlim(1, 97)
        plt.grid(visible=True, color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)

        if id_x == 0:
            legend_handles = [plt.Rectangle((0, 0), 1, 1, color=color) for color in color_groups.values()]
            legend_labels = color_groups.keys()
            plt.legend(legend_handles, legend_labels, bbox_to_anchor=(1.05, 1), loc='upper left')

    return plt

def bootstrap(matrix_file: str, signature_file: str, output_file_exposure_avg: str, output_file_exposure_std: str,
              opportunity_file: str = None, data_type: DataType = DataType.exome):
    """Perform bootstrap analysis of mutational signatures.
    
    Args:
        matrix_file: Path to mutation count matrix file
        signature_file: Path to signature reference file
        output_file_exposure_avg: Path to save average exposure results
        output_file_exposure_std: Path to save exposure standard deviations
        opportunity_file: Optional path to opportunity matrix file
        data_type: Type of data (exome or genome)
    """
    M = read_counts(matrix_file)
    S, index_signature = read_signature(signature_file)
    O = read_opportunity(M, opportunity_file)
    lambd = get_lambda(data_type)
    estimated_exposure, exposure_std = _bootstrap(M, S, O, lambd=lambd)
    np.savetxt(output_file_exposure_avg, estimated_exposure, delimiter='\t')
    np.savetxt(output_file_exposure_std, exposure_std, delimiter='\t')

def download_reference_genome(ref_genome=None, genome_path=None, dest_dir='genomes'):
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
          n_iterations: int = 1000,
          data_type: DataType = DataType.exome):
    """Perform mutational signature refitting analysis.
    
    Args:
        matrix_file: Path to mutation count matrix file
        signature_file: Path to signature reference file
        ref_genome: Reference genome identifier (e.g., GRCh37, GRCh38)
        genome_path: Path to local reference genome file
        n_bootstraps: Number of bootstrap iterations
        opportunity_file: Path to weight matrix file
        numeric_chromosomes: Whether chromosome names are numeric
        genotyped: Whether VCF has genotype information
        output_folder: Directory to save output files
        signature_names: Comma-separated list of signatures to include
        n_iterations: Number of refitting iterations
        data_type: Type of data (exome or genome)
    """
    sig_name = Path(signature_file).stem
    matrix_name = Path(matrix_file).stem
    run_name = f'{matrix_name}_{sig_name}_{n_iterations}'
    logger.info(f'Refitting mutational signatures with {n_iterations} iterations')
    logger.info(f'Run name: {run_name}')
    
    start_time = time.time()
    file_name, file_extension = os.path.splitext(matrix_file)
    
    # Handle VCF input (both .vcf and .vcf.gz)
    if file_extension == '.vcf' or (file_extension == '.gz' and os.path.splitext(file_name)[1] == '.vcf'):
        if not ref_genome and not genome_path:
            raise ValueError("Either ref_genome or genome_path must be provided.")
        genome_path = download_reference_genome(ref_genome=ref_genome, genome_path=genome_path)
        logger.info(f"Reference genome path: {genome_path}")
        # Derive output csv path from input sample name
        sample_name = Path(matrix_file).name
        if sample_name.endswith('.vcf.gz'):
            sample_name = sample_name[:-7]
        elif sample_name.endswith('.vcf'):
            sample_name = sample_name[:-4]
        out_csv = f"{output_folder}/{sample_name}.csv"
        os.makedirs(output_folder, exist_ok=True)
        count_mutation(matrix_file, genome_path, out_csv, numeric_chromosomes, genotyped)
        matrix_file = out_csv
    
    # Read input data
    M = read_counts(matrix_file)

    # Remove rows with all zeros (can occur with genotyped VCFs)
    row_sums = M.sum(axis=1)
    if (row_sums == 0).any():
        n_zero_rows = (row_sums == 0).sum()
        logger.warning(f"Removing {n_zero_rows} empty row(s) with zero mutation counts from the matrix")
        M = M[row_sums > 0]

    index_matrix = M.index.values.tolist()
    S = read_signature(signature_file)
    
    # Filter and validate signatures
    index_matrix = M.index.values.tolist()
    normalized_df = M.div(M.sum(axis=1), axis=0)
    col_means = normalized_df.mean(axis=0)
    count_less_than_01 = (col_means < 0.004).sum()
    threshold = 0.25 * len(col_means)
    
    if count_less_than_01 <= threshold:
        zero_contexts = np.array(M.columns)[(M.values < 0.01).any(axis=0)]
        corr_sigs_mask = (S.loc[:, zero_contexts] >= 0.06).any(axis=1)
        signatures = S.loc[~S.index.isin(corr_sigs_mask[corr_sigs_mask].index)]
        S = signatures
        
    if signature_names is not None:
        requested_sigs = signature_names.split(',')
        # Only keep signatures that are both requested and available after filtering
        available_requested = [sig for sig in requested_sigs if sig in S.index]

        if len(available_requested) < 5:
            missing_sigs = [sig for sig in requested_sigs if sig not in S.index]
            logger.warning(f"Only {len(available_requested)} of the requested signatures are available after correlation filtering.")
            logger.warning(f"Missing signatures: {', '.join(missing_sigs)}")
            logger.warning(f"Available signatures: {', '.join(S.index.tolist())}")
            raise ValueError(f"Only {len(available_requested)} requested signatures are available after filtering. At least 5 are required. Missing: {missing_sigs}")

        if len(available_requested) < len(requested_sigs):
            missing_sigs = [sig for sig in requested_sigs if sig not in S.index]
            logger.warning(f"Some requested signatures were filtered out due to low correlation with the sample: {', '.join(missing_sigs)}")
            logger.info(f"Using {len(available_requested)} available signatures: {', '.join(available_requested)}")

        S = filter_signatures(S, available_requested)
        
    # Prepare data for analysis
    index_signature = S.index.values.tolist()
    desired_order = M.columns
    
    # Check for missing columns
    missing_cols = set(desired_order) - set(S.columns)
    if missing_cols:
        raise ValueError("Error: Following columns are not present in the DataFrame or make sure the file provided is Tab separated:", missing_cols)
        
    # Reorder and convert to numpy arrays
    S = S[desired_order]
    S = S.to_numpy().astype(float)
    M = M.to_numpy().astype(float)
    O = read_opportunity(M, opportunity_file)
    lambd = get_lambda(data_type)
    
    # Handle single sample case
    if (M.ndim == 2 and M.shape[0] == 1) or M.ndim == 1:
        print("The Hyperparameter is", lambd)
        boostrap_M = _bootstrap(M, n_bootstraps)
        #print(boostrap_M)
        logger.info(f"Lambda: {lambd}")
        
        expo_run = _refit(boostrap_M, S, O, lambd=lambd, n_iterations=n_iterations)
        #print(expo_run)
        expo_run = pd.DataFrame(data=expo_run, columns=index_signature)
        mean_columns_E = expo_run.mean(axis=0)
        
        # Adjust threshold based on mutation count
        threshold = 0.06 if any(row.sum() > 600 for row in M) else 0.03
        selected_headers = mean_columns_E[mean_columns_E >= threshold].index.tolist()
        
        # Refit with selected signatures
        S = read_signature(signature_file).loc[selected_headers]
        print(S.shape)
        index_signature = S.index.values.tolist()
        expo_run = _refit(boostrap_M, S, O, lambd=lambd, n_iterations=n_iterations)
        #print(expo_run)
        expo_run = pd.DataFrame(data=expo_run, columns=index_signature)
        
        # Save results
        plot = single_plot(expo_run)
        expo_run_median = expo_run.median()
        expo_run_median.to_csv(f"{output_folder}/StarSign_exposure_median_{run_name}.txt", index=True, sep='\t')
        expo_run.to_csv(f"{output_folder}/StarSign_exposure_Exposure_{run_name}.txt", index=True, header=True, sep='\t')
        
        if plot is not None:
            plot.savefig(f"{output_folder}/StarSign_exposure_Exposure_{run_name}.png", dpi=300)
            plot.close()
        else:
            logger.warning("No plot was generated due to filtering criteria.")
    
    # Handle multiple sample case
    else:
        print("The Hyperparameter is", lambd)
        logger.info(f"Original signature shape: {S.shape}")
        E = _refit(M, S, O, lambd=lambd, n_iterations=n_iterations)
        E = pd.DataFrame(data=E, columns=index_signature, index=index_matrix)
        
        # Filter signatures based on mean exposure
        mean_columns_E = E.mean(axis=0)
        threshold = 0.05 if any(row.sum() > 600 for row in M) else 0.03
        selected_headers = mean_columns_E[mean_columns_E >= threshold].index.tolist()
        
        # Check if any signatures meet the threshold
        if not selected_headers:
            logger.warning(f"No signatures meet the threshold of {threshold}. Using top 5 signatures instead.")
            # Fall back to top 5 signatures by mean exposure
            selected_headers = mean_columns_E.nlargest(5).index.tolist()
        
        # Refit with selected signatures
        S = read_signature(signature_file).loc[selected_headers]
        index_signature = S.index.values.tolist()
        E = _refit(M, S, O, lambd=lambd, n_iterations=n_iterations)
        
        # Calculate and save average exposures
        sum_expo = E.sum(axis=0, keepdims=True) / len(E)
        E = pd.DataFrame(data=E, columns=index_signature, index=index_matrix)
        E.to_csv(f'{output_folder}/{run_name}_threshold.txt', index=index_matrix, header=True, sep='\t')
        
        # Save and plot top signatures
        sum_expo = pd.DataFrame(data=sum_expo, columns=index_signature, index=['Signatures'])
        sum_expo.to_csv(f'{output_folder}/average_{run_name}.txt', columns=index_signature, index=False, sep='\t')
        sum_expo = np.transpose(sum_expo)
        sort_E = sum_expo.sort_values(by=['Signatures'], ascending=False)
        sort_E = sort_E.iloc[:5, 0:]
        index_sort_E = sort_E.index.values.tolist()
        sort_E = pd.DataFrame(data=sort_E, index=index_sort_E)
        
        # Generate plots only if we have data
        if len(sort_E) > 0:
            plot_top_five = cohort_plot(sort_E)
            if plot_top_five is not None:
                plot_top_five.savefig(f"{output_folder}/starsign_top5_signatures_{run_name}.png", dpi=300)
            else:
                logger.warning("Unable to generate top 5 signatures plot.")
        else:
            logger.warning("No signatures available for top 5 plot generation.")
            
        if len(E.columns) > 0:
            plot_variance = cohort_violin(E)
            if plot_variance is not None:
                plot_variance.savefig(f"{output_folder}/starsign_cohort_{run_name}.png", dpi=300)
            else:
                logger.warning("Unable to generate cohort violin plot.")
        else:
            logger.warning("No signatures available for cohort violin plot generation.")
    
    logger.info(f"Analysis completed in {time.time() - start_time:.2f} seconds")

def get_lambda(data_type: DataType) -> float:
    """Get lambda regularization parameter based on data type.
    
    Args:
        data_type: Type of data (exome or genome)
    Returns:
        float: Lambda value for regularization
    """
    return 1000 if data_type == DataType.genome else 0.7


def read_opportunity(M: np.ndarray, opportunity_file: Optional[str] = None) -> np.ndarray:
    """Read or generate opportunity matrix for mutation analysis.
    
    Args:
        M: Mutation count matrix
        opportunity_file: Optional path to opportunity matrix file or predefined type
    Returns:
        np.ndarray: Opportunity matrix
    """
    n_samples = len(M)
    n_mutations = M.shape[1]
    
    if opportunity_file is None:
        return np.ones((n_samples, n_mutations), dtype=float)
        
    if opportunity_file == "human-genome":
        O = np.array([1.14e+08, 6.60e+07, 1.43e+07, 9.12e+07] * 24)  # Simplified array for readability
    elif opportunity_file == "human-exome":
        O = np.array([1940794, 1442408, 514826, 1403756] * 24)  # Simplified array for readability
    else:
        O = pd.read_csv(opportunity_file, sep='\t', header=None).to_numpy().astype(float)
        
    O = np.broadcast_to(O, M.shape)
    O = O / np.amin(O).sum(axis=-1, keepdims=True)
    
    assert O.shape == (n_samples, n_mutations), f'{O.shape} != {(n_samples, n_mutations)}'
    return O

def read_signature(signature_file: str) -> pd.DataFrame:
    """Read signature reference matrix from file.
    
    Args:
        signature_file: Path to signature reference file
    Returns:
        pd.DataFrame: Signature reference matrix
    """
    return pd.read_csv(signature_file, delimiter='\t')

def read_counts(matrix_file: str) -> pd.DataFrame:
    """Read mutation count matrix from file.
    
    Args:
        matrix_file: Path to mutation count matrix file
    Returns:
        pd.DataFrame: Mutation count matrix
    """
    return pd.read_csv(matrix_file, delimiter='\t')

def filter_signatures(S: pd.DataFrame, signature_names: list) -> pd.DataFrame:
    """Filter signature matrix to include only specified signatures.
    
    Args:
        S: Signature reference matrix
        signature_names: List of signature names to include
    Returns:
        pd.DataFrame: Filtered signature matrix
    Raises:
        ValueError: If fewer than 5 signatures are provided
    """
    if len(signature_names) < 5:
        raise ValueError("You must select at least 5 signature names.")
    return S.loc[signature_names]

def get_tri_context_fraction(mut_counts: pd.DataFrame) -> pd.DataFrame:
    """Calculate trinucleotide context fractions from mutation counts.
    
    Args:
        mut_counts: Mutation count matrix
    Returns:
        pd.DataFrame: Trinucleotide context fractions
    """
    total_mutations = mut_counts.sum().sum()
    return mut_counts.div(total_mutations) if total_mutations != 0 else mut_counts

def get_num_cpus() -> int:
    """Get number of available CPU cores.
    
    Returns:
        int: Number of CPU cores
    """
    return multiprocessing.cpu_count()

def denovo(matrix_file: Annotated[str, typer.Argument(help='Tab separated matrix file')],
           n_signatures: Annotated[int, typer.Argument(help='Number of signatures to identify')],
           lambd: Annotated[float, typer.Option(help='Regularization parameter')] = 0.7,
           ref_genome: Annotated[Optional[str], typer.Option(help='Reference genome to use')] = None,
           genome_path: Annotated[Optional[str], typer.Option(help='Local path to reference genome')] = None,
           opportunity_file: Annotated[str, typer.Option(help="Triplet distribution reference")] = None,
           cosmic_file: Annotated[str, typer.Option(help='COSMIC signature reference')] = None,
           numeric_chromosomes: Annotated[bool, typer.Option(help="Numeric chromosome names")] = False,
           genotyped: Annotated[bool, typer.Option(help="VCF has genotype information")] = True,
           max_em_iterations: int = 100,
           max_gd_iterations: int = 50,
           output_folder: str = 'output/'):
    """Perform de novo mutational signature analysis.
    
    Args:
        matrix_file: Path to mutation count matrix
        n_signatures: Number of signatures to identify
        lambd: Regularization parameter
        ref_genome: Reference genome identifier
        genome_path: Path to reference genome file
        opportunity_file: Path to opportunity matrix
        cosmic_file: Path to COSMIC signature file
        numeric_chromosomes: Whether chromosome names are numeric
        genotyped: Whether VCF has genotype information
        max_em_iterations: Maximum EM algorithm iterations
        max_gd_iterations: Maximum gradient descent iterations
        output_folder: Directory to save outputs
    """
    matrix_name = Path(matrix_file).stem
    run_name = matrix_name
    logger.info(f'Starting de novo analysis for {run_name}')
    start_time = time.time()

    # Handle VCF input (both .vcf and .vcf.gz)
    if matrix_file.endswith('.vcf') or matrix_file.endswith('.vcf.gz'):
        if not ref_genome and not genome_path:
            raise ValueError("Either ref_genome or genome_path must be provided.")
        genome_path = download_reference_genome(ref_genome=ref_genome, genome_path=genome_path)
        logger.info(f"Reference genome path: {genome_path}")
        # Derive output csv path from input sample name
        sample_name = Path(matrix_file).name
        if sample_name.endswith('.vcf.gz'):
            sample_name = sample_name[:-7]
        elif sample_name.endswith('.vcf'):
            sample_name = sample_name[:-4]
        out_csv = f"{output_folder}/{sample_name}.csv"
        os.makedirs(output_folder, exist_ok=True)
        count_mutation(matrix_file, genome_path, out_csv, numeric_chromosomes, genotyped)
        matrix_file = out_csv

    # Read and prepare data
    M = read_counts(matrix_file)

    # Remove rows with all zeros (can occur with genotyped VCFs)
    row_sums = M.sum(axis=1)
    if (row_sums == 0).any():
        n_zero_rows = (row_sums == 0).sum()
        logger.warning(f"Removing {n_zero_rows} empty row(s) with zero mutation counts from the matrix")
        M = M[row_sums > 0]

    index_matrix = M.index.values.tolist()
    desired_order = M.columns
    O = read_opportunity(M, opportunity_file)
    M = M.to_numpy().astype(float)

    # Perform de novo analysis
    E, S = _denovo(M, n_signatures, lambd, O, em_steps=max_em_iterations, gd_steps=max_gd_iterations)

    # Compare with COSMIC signatures if provided
    if cosmic_file is not None:
        cosmic = pd.read_csv(cosmic_file, delimiter='\t')
        missing_cols = set(desired_order) - set(cosmic.columns)
        if missing_cols:
            raise ValueError("Error: Following columns are not present in the DataFrame:", missing_cols)
        cosmic = cosmic[desired_order]
        cos_similarity = cos_sim_matrix(S, cosmic)[0]
        cos_similarity.to_csv(f"{output_folder}/StarSign_{run_name}_cosine_similarity.txt", sep="\t")

    # Prepare and save results
    S = np.transpose(S)
    Sig = ['De novo ' + letter for letter in string.ascii_uppercase[:S.shape[1]]]
    mutation_labels = ['A[C>A]A', 'A[C>A]C', 'A[C>A]G', 'A[C>A]T'] * 24  # Simplified for readability
    
    S = pd.DataFrame(S, columns=Sig, index=mutation_labels)
    S.to_csv(f"{output_folder}/StarSign_{run_name}_Denovo_signature.txt", sep='\t', index=mutation_labels)
    
    # Generate and save plots
    deno_figure = plot_profile(S.T)
    deno_figure.savefig(f"{output_folder}/StarSign_{run_name}_profile.png", dpi=600)
    
    E = pd.DataFrame(E, columns=Sig, index=None)
    E.to_csv(f"{output_folder}/StarSign_{run_name}_Denovo_exposures.txt", sep='\t', index=None)
    
    logger.info(f"Analysis completed in {time.time() - start_time:.2f} seconds")

def main():
    """Main entry point for the CLI application."""
    app = typer.Typer()
    app.command()(refit)
    app.command()(denovo)
    app.command()(count_mutation)
    app()

if __name__ == "__main__":
    main()
