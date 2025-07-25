========
StarSignDNA: Signature tracing for accurate representation of mutational processes
========

A Robust Tool for Mutational Signature Analysis
---------------------------------------------

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
    :target: https://opensource.org/licenses/MIT
    :alt: MIT License

.. image:: https://img.shields.io/badge/release-v1.0.0-brightgreen.svg
    :target: https://pypi.org/project/starsigndna/
    :alt: Release Version

StarSignDNA is a powerful algorithm for mutational signature analysis that provides both efficient refitting of known signatures and de novo signature extraction. The tool excels at:

* Deciphering well-differentiated signatures linked to known mutagenic mechanisms
* Demonstrating strong associations with patient clinical features
* Providing robust statistical analysis and visualization
* Handling both single-sample and cohort analyses
* Supporting various input formats including VCF files

Quick Links
----------
* **License**: `MIT license <https://opensource.org/licenses/MIT>`_
* **Paper**: `Preprint <https://www.biorxiv.org/content/10.1101/2024.06.29.601345v1>`_
* **Source Code**: https://github.com/uio-bmi/StarSignDNA

Installation
-----------

You can install StarSignDNA in one of two ways:

1. Via PyPI (Recommended)::

    pip install starsigndna

2. From Source::

    git clone https://github.com/uio-bmi/StarSignDNA
    cd StarSignDNA
    pip install -e .

Key Features
-----------

* **Signature Refitting**: Accurate refitting of known mutational signatures to new samples
* **De Novo Discovery**: Novel signature extraction from mutation catalogs
* **Flexible Input**: Support for both VCF files and mutation count matrices
* **Statistical Robustness**: Built-in bootstrapping and cross-validation
* **Visualization**: Comprehensive plotting functions for signature analysis
* **Performance**: Optimized algorithms for both small and large datasets
* **Reference Genome Support**: Built-in handling of various genome builds
* **Enhanced CLI**: Improved signature parsing supporting both comma and space-separated formats

Basic Usage
----------

Get help on available commands::

    starsigndna --help

Available Commands:
~~~~~~~~~~~~~~~~~

* **count-mutation**: Count mutation types in VCF files
* **denovo**: Perform de novo signature discovery
* **refit**: Refit known signatures to new samples

Signature Refitting
------------------

The refitting algorithm matches mutation patterns against known COSMIC signatures.

Basic Usage::

    starsigndna refit <matrix_file> <signature_file> [OPTIONS]

Example with Specific Signatures::

    starsigndna refit example_data/M_catalogue.txt example_data/COSMICv34.txt \
        --output-folder /test_result \
        --signature-names "SBS40c,SBS2,SBS94"

Example with Space-Separated Signatures::

    starsigndna refit example_data/M_catalogue.txt example_data/COSMICv34.txt \
        --output-folder /test_result \
        --signature-names "SBS40c SBS2 SBS94"

Example with VCF Input::

    starsigndna refit example_data/tcga_coad_single.vcf example_data/sig_cosmic_v3_2019.txt \
        --output-folder /output \
        --signature-names "SBS40c,SBS2,SBS94" \
        --ref-genome GRCh37

Key Options:
~~~~~~~~~~~

* **--ref_genome**: Reference genome for VCF processing
* **--n_bootstraps**: Number of bootstrap iterations (default: 200)
* **--opportunity_file**: Custom mutation opportunity matrix
* **--signature_names**: Specific signatures to consider (minimum 5 signatures required)
* **--n_iterations**: Maximum optimization iterations (default: 1000)

**Signature Names Format**: The `--signature-names` parameter accepts both comma-separated and space-separated formats:
* Comma-separated: `"SBS1,SBS3,SBS5,SBS6,SBS8"`
* Space-separated: `"SBS1 SBS3 SBS5 SBS6 SBS8"`

Expected Output:
~~~~~~~~~~~~~~~

The refit command generates several output files in the specified output folder:

**For single sample analysis:**
* **StarSign_refit_exposure_median_{run_name}.txt**: Median exposure values across bootstrap iterations
* **StarSign_refit_exposure_Exposure_{run_name}.txt**: Full exposure matrix from bootstrap analysis
* **StarSign_refit_exposure_Exposure_{run_name}.png**: Violin plot of exposure distributions

**For cohort analysis:**
* **refit_{run_name}_threshold.txt**: Exposure matrix after signature filtering
* **average_refit_{run_name}.txt**: Average exposure values across samples
* **starsign_refit_top5_signatures_{run_name}.png**: Bar plot of top 5 signatures by average exposure
* **starsign_refit_cohort_{run_name}.png**: Violin plot showing exposure distributions across cohort

**For VCF input:**
* **matrix.csv**: Generated mutation count matrix from VCF file

De Novo Signature Discovery
-------------------------

The de novo algorithm extracts novel signatures from mutation patterns.

Basic Usage::

    starsigndna denovo <matrix_file> <n_signatures> [OPTIONS]

Example with Optimal Parameters::

    starsigndna denovo snakemake/results/data/M_catalogue.txt 4 0.1 \
        --cosmic-file example_data/COSMICv34.txt \
        --output-folder /test_result

Parameter Optimization
~~~~~~~~~~~~~~~~~~~~

1. Configure Grid Search Parameters::

    cd snakemake
    vi Snakefile

    # Example configuration:
    ks = list(range(2, 10))  # Number of signatures
    lambdas = [0, 0.01, 0.05, 0.1, 0.2]  # Regularization values

2. Run Grid Search::

    snakemake -j <num_cpu>

3. Find Optimal Parameters::

    sort -k3n,3 results/data/all.csv

Key Options:
~~~~~~~~~~~

* **--lambd**: Regularization parameter (default: 0.7)
* **--opportunity-file**: Custom mutation opportunity matrix
* **--cosmic-file**: Reference signatures for comparison
* **--max-em-iterations**: Maximum EM iterations (default: 100)
* **--max-gd-iterations**: Maximum gradient descent steps (default: 50)

Expected Output:
~~~~~~~~~~~~~~~

The denovo command generates several output files in the specified output folder:

* **StarSign_denovo_{run_name}_signature.txt**: Extracted mutational signatures matrix (signatures × mutation types)
* **StarSign_denovo_{run_name}_exposures.txt**: Signature exposures for each sample (samples × signatures)
* **StarSign_denovo_{run_name}_profile.png**: Visualization of extracted signatures
* **StarSign_denovo{run_name}_cosine_similarity.txt**: Similarity scores with known COSMIC signatures (if --cosmic-file provided)

**For VCF input:**
* **matrix.csv**: Generated mutation count matrix from VCF file

Advanced Features
---------------

* **Opportunity Matrices**: Support for custom mutation opportunity matrices:
  - Built-in human genome/exome distributions
  - Custom tissue-specific distributions
  - Patient-specific normal tissue references

* **Input Flexibility**: 
  - VCF files (single or multi-sample)
  - Pre-computed mutation matrices
  - Various chromosome naming conventions

* **Output Customization**:
  - Detailed signature profiles
  - Exposure matrices
  - Visualization plots
  - Statistical metrics

* **Enhanced Error Handling**:
  - Robust signature filtering with fallback mechanisms
  - Graceful handling of empty datasets
  - Improved plotting error recovery

Recent Improvements
------------------

* **Enhanced CLI**: Improved signature parsing supporting both comma and space-separated formats
* **Better Error Handling**: Robust signature filtering with automatic fallback to original signature set
* **Improved Plotting**: Enhanced visualization functions with better error recovery
* **Code Quality**: Comprehensive comments and documentation added to all scripts
* **Snakemake Integration**: Enhanced workflow scripts with better reproducibility and error handling

Contributing
-----------

We welcome contributions! Please feel free to submit a Pull Request.

Contact
-------

* **Maintainer**: Christian D. Bope
* **Email**: christianbope@gmail.com
* **Institution**: University of Oslo

Citation
--------

If you use StarSignDNA in your research, please cite our paper:
[Citation details to be added after publication]


