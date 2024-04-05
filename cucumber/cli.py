"""Console script for cucumber."""
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


def single_plot(df):
    header = df.columns.tolist()

# 2. Compute the median per column and output as a vector with the row name 'signatures'
    sorted_df = df.sort_values(by=list(df.columns))
#    print(sorted_df)
    median = sorted_df.median()
    median.name = 'signatures'

# 3. Compute the 2.5 percentile of the dataframe and output the vector with the row name 'perc2.5'
#perc_2_5 = df.quantile(0.025)
    perc_2_5 = df.quantile(0.25)
    perc_2_5.name = 'perc2.5'

# 4. Compute the 97.5 percentile of the dataframe and output the vector with the row name 'perc97.5'
#perc_97_5 = df.quantile(0.975)
    perc_97_5 = df.quantile(0.75)
    perc_97_5.name = 'perc97.5'

# 5. Append the 3 vectors in a dataframe with the header and the following row names ['Signature', 'E_25', 'E_95']
    output_df = pd.concat([median, perc_2_5, perc_97_5], axis=1).transpose()
    output_df.index = ['Signature', 'E_25', 'E_75']
    output_df.columns = header
    data = output_df.transpose()
    filtered_data = data[data['Signature'] >= 0.01]
#    print(filtered_data)
    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(20, 8))
    for i, row in filtered_data.iterrows():
        lower_error = row['E_25']
        upper_error = row['E_75']
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


def refit(matrix_file: Annotated[str, typer.Argument(help='Tab separated matrix file')],
          signature_file: Annotated[str, typer.Argument(help='Comma separated matrix file')],
          opportunity_file: str = None, ref_genome: str = None, n_bootstraps: int = 200,
          numeric_chromosomes: bool = None, genotyped: bool = None, cancer_type: str = None,
          output_folder: str = 'output/'):
    #    numeric_chromosomes: Annotated[bool, typer.Argument(help="True if chromosome names in vcf are '1', '2', '3'. False if 'chr1', 'chr2', 'chr3'")] = True,
    #   genotyped: Annotated[bool, typer.Argument(help="True if the VCF file has genotype information for many samples")] = False, output_folder: str = 'output/',
    #   cancer_type: Annotated[str,typer.Argument(help="Cancer type abbreviation, eg.: bcla, brca, chol, gbm, lgg, cesc, coad, esca, uvm, hnsc, kich, kirp, kirc, lihc, luad, lusc, dlbc, laml, ov, paad, prad, sarc, skcm, stad, thca, ucec")] = None):
    '''
    Mutational SIgnatures Refit Parameters \n
    n_bootstraps: \n
    matrix_file: str \n
    signature_file: str \n
    output_file_exposure: str \n
    opportunity_file: str \n
    numeric_chromosomes,bool: True if chromosome names in vcf are '1', '2', '3'. False if 'chr1', 'chr2', 'chr3' \n
    genotyped,bool:True if the VCF file has genotype information for many samples \n
    cancer_type, str: bcla, brca, chol, gbm, lgg, cesc, coad, esca, uvm, hnsc, kich, kirp, kirc, lihc, luad, lusc, dlbc, laml, ov, paad, prad, sarc, skcm, stad, thca, ucec

    '''
    start_time = time.time()
    # reference genome
    #ref_genome = '/Users/bope/Documents/MutSig/scientafellow/packages/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
    file_name, file_extension = os.path.splitext(matrix_file)
    if file_extension == '.vcf':
        assert ref_genome is not None, 'Please provide a reference genome along with the vcf file'
        count_mutation(matrix_file, ref_genome, f'{output_folder}/matrix.csv', numeric_chromosomes, genotyped)
        matrix_file = f'{output_folder}/matrix.csv'
    M, index_matrix = read_counts(matrix_file)
    #print("MMMM",M.columns)
    S, index_signature = read_signature(signature_file)
    desired_order = M.columns

# Check for missing columns
    missing_cols = set(desired_order) - set(S.columns)  # Set difference for missing columns
    if len(missing_cols) > 0:
        raise ValueError("Error: Following columns are not present in the DataFrame:", missing_cols)
# Reorder columns
    S = S[desired_order]
    S = S.to_numpy().astype(float)
    M = M.to_numpy().astype(float)
    if cancer_type is not None:
        true_order = 'Type	SBS1	SBS2	SBS3	SBS4	SBS5	SBS6	SBS7a	SBS7b	SBS7c	SBS7d	SBS8	SBS9	SBS10a	SBS10b	SBS10c	SBS10d	SBS11	SBS12	SBS13	SBS14	SBS15	SBS16	SBS17a	SBS17b	SBS18	SBS19	SBS20	SBS21	SBS22a	SBS22b	SBS23	SBS24	SBS25	SBS26	SBS27	SBS28	SBS29	SBS30	SBS31	SBS32	SBS33	SBS34	SBS35	SBS36	SBS37	SBS38	SBS39	SBS40a	SBS40b	SBS40c	SBS41	SBS42	SBS43	SBS44	SBS45	SBS46	SBS47	SBS48	SBS49	SBS50	SBS51	SBS52	SBS53	SBS54	SBS55	SBS56	SBS57	SBS58	SBS59	SBS60	SBS84	SBS85	SBS86	SBS87	SBS88	SBS89	SBS90	SBS91	SBS92	SBS93	SBS94	SBS95	SBS96	SBS97	SBS98	SBS99'.split()[
                     1:]
        #  assert index_signature.tolist() != true_order, (f'The order of the signatures in the signature file is not the same as cosmic 3.4. Cannot do automatic selection of {cancer_type} signatures', index_signature.tolist(), true_order)
        assert index_signature == true_order, (
        f'The order of the signatures in the signature file is not the same as cosmic 3.4. You can download the file at https://cancer.sanger.ac.uk/cosmic/download/cosmic. Cannot do automatic selection of {cancer_type} signatures',
        index_signature, true_order)
        if cancer_type == 'bcla':
            index = [0, 1, 3, 4, 18]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'brca':
            index = [0, 1, 2, 4, 18, 24, 37, 41, 50]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'chol':
            index = [0, 1, 4, 18, 24, 47, 48, 49, 67]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'gbm':
            index = [0, 4, 37, 47, 48, 49]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'lgg':
            index = [0, 4]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'cesc':
            index = [0, 1, 4, 18]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'coad':
            index = [0, 4, 20, 47, 48, 49, 53]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'esca':
            index = [0, 2, 4, 18, 22, 23, 24, 47, 48, 49]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'uvm':
            index = [0, 4, 47, 48, 49, 60, 61]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'hnsc':
            index = [0, 1, 3, 4, 18, 47, 48, 49, 54]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'kich':
            index = [0, 1, 18, 36, 47, 48, 49]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'kirp':
            index = [0, 1, 4, 18, 54, 58]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'kirc':
            index = [0, 4, 28, 29, 47, 48, 49, 50]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'lihc':
            index = [0, 2, 3, 4, 17, 21, 24, 28, 29, 36, 47, 48, 49]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'luad':
            index = [0, 1, 3, 4, 18, 47, 48, 54]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'lusc':
            index = [0, 1, 3, 4, 18, 54]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'dlbc':
            index = [1, 2, 4, 5, 11, 18, 22, 23, 41, 43, 44, 47, 48, 49, 65]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'dlbc':
            index = [0, 1, 2, 4, 5, 11, 18, 22, 23, 41, 43, 44, 47, 48, 49]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'laml':
            index = [0, 4, 52]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'ov':
            index = [0, 1, 2, 4, 18, 24, 47, 48, 49]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'paad':
            index = [0, 1, 2, 4, 18, 22, 23, 24, 47, 48, 49]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'prad':
            index = [0, 2, 4, 24, 47, 48, 49, 67]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'sarc':
            index = [0, 4, 6, 7, 8, 9, 24, 68]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'skcm':
            index = [0, 4, 6, 7, 8, 9, 45, 47, 48, 49, 54, 58]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'stad':
            index = [0, 1, 2, 18, 20, 22, 23, 24, 26, 47, 48, 49]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'thca':
            index = [0, 1, 4, 18, 47, 48, 49, 52, 54, 58, 62]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        elif cancer_type == 'ucec':
            index = [0, 1, 4, 12, 13, 14, 15, 18, 19, 20, 35, 53]
            S = S[index]
            index_signature = [index_signature[i] for i in index]
        else:
            raise ValueError(
                f'Unknown cancer type {cancer_type}. Valid cancer types are: bcla, brca, chol, gbm, lgg, cesc, coad, esca, uvm, hnsc, kich, kirp, kirc, lihc, luad, lusc, dlbc, laml, ov, paad, prad, sarc, skcm, stad, thca, ucec')
    O = read_opportunity(M, opportunity_file)
    lambd = 0.7
    print(M.shape)
    if (M.ndim == 2 and M.shape[0] == 1) or M.ndim == 1:
        print("HERRRRRRE")
   #     E = _refit(M, S, O, lambd=lambd)[0]
        boostrap_M = _bootstrap(M,n_bootstraps)
    #    print(E.shape)
        expo_run = _refit(boostrap_M, S, O, lambd=lambd)
 #       expo_run = _bootstrap(M, S, O, n_bootstraps, lambd=lambd)
 ##       E = [E, E_std]
    #    E = [E,E_25,E_95]
    #    E = pd.DataFrame(data=E, columns=index_signature, index=['Signature', 'E_25', 'E_95'])
        expo_run = pd.DataFrame(data=expo_run, columns=index_signature)
      #  print(expo_run)
        plot = single_plot(expo_run)
        #expo_run = np.array(expo_run)
        np.savetxt(f'{output_folder}/boostrap_catalogue_test.txt', np.array(boostrap_M))
    #    np.savetxt(f'{output_folder}/Exposure_3avril.txt', np.array(E))
 #       np.savetxt(f'{output_folder}/Exposure_run_3avril_1.txt', np.array(expo_run))
        expo_run.to_csv(f"{output_folder}/StarSign_exposure_Exposure_test.txt", index=True, header=True, sep='\t')
        plot.savefig(f"{output_folder}/StarSign_exposure_Exposure_test.png", dpi=600)
    else:
        print("NOOOOOOOOO")
        E = _refit(M, S, O, lambd=lambd)
        # print(O)
        sum_expo = E.sum(axis=0, keepdims=True) / len(E)
        sum_expo_t = np.transpose(sum_expo)
        E = pd.DataFrame(data=E, columns=index_signature, index=index_matrix)
        E.to_csv(f'{output_folder}/sim5_cosmic67_binomial_common_500_l07_l0_19mars_5000.txt', index=index_matrix, header=True, sep='\t')
        sum_expo = pd.DataFrame(data=sum_expo, columns=index_signature, index=['Signatures'])
        sum_expo = np.transpose(sum_expo)
        plot_summary = cohort_plot(sum_expo)
        plot_summary.savefig(f"{output_folder}/sim5_cosmic67_binomial_common_500_l07_l0_19mars_5000.png", dpi=600)
        sort_E = sum_expo.sort_values(by=['Signatures'], ascending=False)
        sort_E = sort_E.iloc[:5, 0:]
        plot_top_five = cohort_plot(sort_E)
        plot_top_five.savefig(f"{output_folder}/sim5_cosmic67_binomial_common_500_l07_l0_19mars_5000.png", dpi=600)
        plot_variance = cohort_violin(E)
        plot_variance.savefig(f"{output_folder}/sim5_cosmic67_binomial_common_500_l07_l0_19mars_5000.png", dpi=600)
        # sum_expo.to_csv(f'{output_folder}/average_exposure_cohort.txt', index=index_signature,
        #                 header=True,
        #                 sep='\t')
        np.savetxt(f'{output_folder}/average_sim5_cosmic67_binomial_common_500_l07_l0_19mars_5000.txt', np.array(sum_expo_t))
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

    # The exposure is normalize to have the same proportion to the catalogue matrix
    #     normalized_vector1 = O / np.linalg.norm(O)
    #     min_value_vector2 = np.min(M)
    #     max_value_vector2 = np.max(M)
    #     O = normalized_vector1 * (max_value_vector2 - min_value_vector2) + min_value_vector2

    assert O.shape == (n_samples, n_mutations), f'{O.shape} != {(n_samples, n_mutations)}'
    return O


def read_signature(signature_file):
    S = pd.read_csv(signature_file, delimiter='\t')
    #S = S.T
    index_signature = S.index.values.tolist()
    #S = S.to_numpy().astype(float)
    # S = S[0:10,0:96]
    return S, index_signature


def read_counts(matrix_file):
    M = pd.read_csv(matrix_file, delimiter='\t')
  #  M = M.T
    index_matrix = M.index.values.tolist()
    #M = M.to_numpy().astype(float)
    return M, index_matrix


def get_tri_context_fraction(mut_counts):
    total_mutations = mut_counts.sum().sum()
    trinucleotide_fractions = mut_counts.div(total_mutations) if total_mutations != 0 else mut_counts
    return trinucleotide_fractions


def get_num_cpus():
    return multiprocessing.cpu_count()


def denovo(matrix_file: Annotated[str, typer.Argument(help='Tab separated matrix file')],
           n_signatures: int, lambd: float = 0.7, #Annotated[float, typer.Argument(help='Regularization parameter')] = 0.7,
  ###         lambd: Annotated[float, typer.Argument(help='Regularization parameter')] = 0.7,
           opportunity_file: str = None,
           cosmic_file: str= None, #Annotated[str, typer.Argument(help='Comma separated cosmic file')] = None,
    ##       cosmic_file: Annotated[str, typer.Option(help='Comma separated cosmic file')] = None,
           numeric_chromosomes: Annotated[bool, typer.Argument(help="True if chromosome names in vcf are '1', '2', '3'. False if 'chr1', 'chr2', 'chr3'")] = False,
           genotyped: Annotated[bool, typer.Argument(help="True if the VCF file has genotype information for many samples")] = True,
           max_em_iterations: int = 100,
           max_gd_iterations: int = 50,
   #        numeric_chromosomes: Annotated[bool, typer.Argument(help="True if chromosome names in vcf are '1', '2', '3'. False if 'chr1', 'chr2', 'chr3'")] = False,
   #        genotyped: Annotated[bool, typer.Argument(help="True if the VCF file has genotype information for many samples")] = True,
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
    start_time = time.time()
    if file_extension == '.vcf':
        assert ref_genome is not None, 'Please provide a reference genome along with the vcf file'
        count_mutation(matrix_file, ref_genome, f'{output_folder}/matrix.csv', numeric_chromosomes, genotyped)
        matrix_file = f'{output_folder}/matrix.csv'
    M, index_matrix = read_counts(matrix_file)
    #print(M)
    n_samples = len(M)
    n_signatures = n_signatures
    lambd = lambd
    print('Lambda',lambd)
    n_mutations = M.shape[1]
    O = read_opportunity(M, opportunity_file)
    # O = O / np.amax(O).sum(axis=-1, keepdims=True)
    # assert O.shape == (n_samples, n_mutations), f'{O.shape} != {(n_samples, n_mutations)}'
    desired_order = M.columns
    #print(desired_order)
    M = M.to_numpy().astype(float)
    E, S = _denovo(M, n_signatures, lambd, O, em_steps=max_em_iterations, gd_steps=max_gd_iterations)
    if cosmic_file is not None:
        cosmic = pd.read_csv(cosmic_file, delimiter='\t')
        #desired_order = M.columns

# Check for missing columns
        missing_cols = set(desired_order) - set(cosmic.columns)  # Set difference for missing columns
        if len(missing_cols) > 0:
            raise ValueError("Error: Following columns are not present in the DataFrame:", missing_cols)
# Reorder columns
        cosmic = cosmic[desired_order]
        cos_similarity = cos_sim_matrix(S, cosmic)[0]
        cos_similarity.to_csv(f"{output_folder}/StarSign_Cosine_similarity_denovo.txt", sep="\t")
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
    S.to_csv(f"{output_folder}/StarSign_Denovo_signature.txt", sep='\t', index=label)
    # print(S)
    deno_figure = plot_profile(S)
    deno_figure.savefig(f"{output_folder}/StarSign_denovo.png", dpi=600)
    E = pd.DataFrame(E, columns=Sig, index=None)
    # np.savetxt(output_file_exposure, np.array(E))
    E.to_csv(f"{output_folder}/StarSign_Denovo_exposures.txt", sep='\t', index=None)
    # np.savetxt(output_file_signature, np.array(S))
    print("--- %s seconds ---" % (time.time() - start_time))


def main():
    app = typer.Typer()
    app.command()(refit)
    app.command()(denovo)
    app.command()(count_mutation)
    app()


if __name__ == "__main__":
    main()
