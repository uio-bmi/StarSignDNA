========
StarSignDNA
========
StarSignDNA algorithm for mutational signature analysis which offers efficient refitting and de novo
mutational signature extraction. StarSignDNA is capable of deciphering well-differentiated signatures linked to known
mutagenic mechanisms and demonstrates strong associations with patient clinical features. The package offers a user-
friendly interface and data visualization routines.

Mutation signature analysis package
-----------------------------------

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
    :target: https://opensource.org/licenses/MIT
    :alt: MIT License

* Free software: `MIT license <https://opensource.org/licenses/MIT>`_
* Documentation: Link to be provided



Paper
--------
Link to be provided

Features
--------

Stable release
--------------

To install StarSign, run this command in your terminal::

    1. Download StarSign from https://github.com/uio-bmi/StarSignDNA
    2. Unzip StarSignDNA-master.zip
    3. cd StarSigndna-master/
    4. pip install -e .

Alternatively, you can install the package directly from PyPi::

    pip install starsigndna
Getting started
---------------

To obtain help::

    starsigndna --help

Usage: starsigndna [OPTIONS] COMMAND [ARGS]...


Options
-------

- **--install-completion**: Install completion for the current shell.
- **--show-completion**: Show completion for the current shell, to copy it or customize the installation.
- **--help**: Show this message and exit.

Commands
--------

- **count-mutation**: Count mutation types in VCF file.
- **denovo**: Performs denovo Mutational Signatures analysis.
- **refit**: Mutational Signatures Refit Parameters

Running mutational signature refit algorithm
---------------------------------------------

The refitting algorithm takes as input a mutational catalog and a COSMIC mutational signature file. The user can also specify signatures to be considered instead of using the full COSMIC matrix or a subset matrix::

    starsigndna refit --help

StarSignDNA Refit Parameters
----------------------------

**Arguments**

- **matrix_file** (TEXT): Tab separated matrix file [default: None] [required]
- **signature_file** (TEXT): Tab separated matrix file [default: None] [required]

**Options**

- **--ref_genome** (TEXT): Path to the reference genome [default: None]
- **--n_bootstraps** (INTEGER): Number of bootstraps [default: 200]
- **--opportunity_file** (TEXT): Path to the opportunity file [default: None]
- **--numeric_chromosomes**: If set, chromosome names are numeric [default: no-numeric-chromosomes]
- **--no_numeric_chromosomes**: If set, chromosome names are not numeric [default: no-numeric-chromosomes]
- **--genotyped**: If set, VCF file has genotype information for many samples [default: genotyped]
- **--no_genotyped**: If set, VCF file does not have genotype information for many samples [default: genotyped]
- **--output_folder** (TEXT): Path to the output folder [default: output/]
- **--signature_names** (TEXT): Comma separated list of signature names [default: None]
- **--n_iterations** (INTEGER): Number of iterations [default: 1000]
- **--help**: Show this message and exit

Running StarSignDNA refit::

    starsigndna refit example_data/M_catalogue.txt example_data/COSMICv34.txt --output-folder /test_result --signature-names SBS40c,SBS2,SBS94
    starsigndna refit example_data/tcga_coad_single.vcf example_data/sig_cosmic_v3_2019.txt --output-folder /output --signature-names SBS40c,SBS2,SBS94 --ref-genome

When the --signature-names option is used, the default number of signature is 3, but we recommend minimum of 5 signatures.
The test data is provided in the example_data folder. To convert *.vcf to a matrix, the user must provide the path to the reference genome using the option --ref-genome.

The user can also provide the distribution of triplets in a reference genome/exome or normal tissue in the same patient (Opportunity matrix) using the option --opportunity-file human-genome/human-exome.


Running mutational signature de novo algorithm
-----------------------------------------------

The de novo algorithm takes as input a mutational catalog and infers the exposure matrix and mutational signature matrix. The COSMIC mutational signature file is provided to compute the cosine similarity::

    starsigndna denovo --help

Performs denovo Mutational Signatures analysis
===============================================

**Arguments**

- **matrix_file** (TEXT): Tab separated matrix file [default: None] [required]
- **n_signatures** (INTEGER): Number of signatures to identify [default: None] [required]

**Options**

- **--lambd** (FLOAT): Regularization parameter [default: 0.7]
- **--opportunity-file** (TEXT): The distribution of triplets in a reference 'human-genome' or 'human-exome' or normal tissue [default: None]
- **--cosmic-file** (TEXT): Tab separated cosmic file [default: None]
- **--numeric-chromosomes**: If set, chromosome names are numeric [default: no-numeric-chromosomes]
- **--no-numeric-chromosomes**: If set, chromosome names are not numeric [default: no-numeric-chromosomes]
- **--genotyped**: If set, VCF file has genotype information for many samples [default: genotyped]
- **--no-genotyped**: If set, VCF file does not have genotype information for many samples [default: genotyped]
- **--max-em-iterations** (INTEGER): Maximum EM iterations [default: 100]
- **--max-gd-iterations** (INTEGER): Maximum GD iterations [default: 50]
- **--file-extension** (TEXT): File extension [default: None]
- **--ref-genome** (TEXT): Path to the reference genome [default: None]
- **--output-folder** (TEXT): Path to the output folder [default: output/]
- **--help**: Show this message and exit

Step 1: Grid Search: The grid uses cross-validation to find the optimal pairwise (k and λ) by going to the snakemake folder and opening the running file (Snakefile) to check all the paths and input files::

    cd snakemake
    vi Snakefile

Step 2: In the Snakefile, provide the range of the number of signatures k and λ for the grid search to determine the optimal k and λ::

    localrules: all
    ks = list(range(2, 10)): default range of the number of signatures
    lambdas = [0, 0.01, 0.05, 0.1, 0.2]: default range of λ

Input mutational catalogue needs to be provided in the dataset folder::

    rule test_train_split:
        input: "results/{dataset}/M_catalogue.txt"

Running the grid search::

    snakemake -j num_cpu

To check manually the optimal k and λ from the output::

    sort -k3n,3 results/data/all.csv

Run denovo using optimal k=4 and λ=0.1::

    starsigndna denovo snakemake/results/data/M_catalogue.txt 4 0.1 --cosmic-file example_data/COSMICv34.txt --output-folder /test_result

Contact
-------

Maintainer Name - chrisbop@uio.no or christianbope@gmail.com


