========
Cucumber
========


.. image:: https://img.shields.io/pypi/v/cucumber.svg
        :target: https://pypi.python.org/pypi/cucumber

.. image:: https://img.shields.io/travis/knutdrand/cucumber.svg
        :target: https://travis-ci.com/knutdrand/cucumber

.. image:: https://readthedocs.org/projects/cucumber/badge/?version=latest
        :target: https://cucumber.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status




Mutation signature analysis package
-----------------------------------


* Free software: MIT license
* Documentation: https://cucumber.readthedocs.io.


Features
--------

Getting started
---------------

Install cucumber by running:
::
   % pip install cucumber
To obtain help :
::
   % cucumber --help
   $ Usage: cucumber [OPTIONS] COMMAND [ARGS]...
Commands:
::
  count-mutation  Count mutation types in a VCF file
  denovo Parameters ---matrix_file: str n_signatures 
  refit  Parameters ---numeric_chromosomes n_bootstraps
  
Running mutational signature refit algorithm:
-----------------------------------------------
The refitting algorithm takes as input a mutational catalog and cosmic mutational signature file
::
  %cucumber refit --help
 
Running cucumber refit
::
  %cucumber refit example_data/skin20.txt example_data/sig_cosmic_v3_2019.txt output/expo.txt output/expo_avg.txt
  %cucumber refit example_data/tcga_coad_single.vcf example_data/sig_cosmic_v3_2019.txt output/expo.txt output/expo_avg.txt
The test data is provided in example_data folder
<<<<<<< HEAD

:: output files for a single sample
::
   $output_file_exposure: exposure matrix with std_dev 
   $exposures_single_dotplot.png: exposure matrix plot with std_dev
The standard deviation is computed using a default of 100 bootstraps. 

.. image:: output/exposures_single_dotplot_skin.png
   :width: 600

:: output files for a cohort
::
   $exposures_cohort_variance: a plot showing the variance of each sample and the mean exposures
   $output_file_exposure: a cohort exposures matrix
   $exposures_cohort_top_5: a plot showing the top 5 exposures
   $ exposures_cohort_dotplot: a plot showing a dotplot of the exposure matrix

.. image:: output/exposures_cohort_variance.png
   :width: 600

.. image:: output/exposures_cohort_top_5.png
   :width: 600
Running mutational signature de novo algorithm:
-----------------------------------------------
The de novo algorithm takes as input a mutational catalog and inferred the exposure matrix and mutational signature matrix. The cosmic mutational signature file is provided to compute the cosine similarity.  
::
  % cucumber denovo --help
:: Step 1: Grid Search: The grid uses cross-validation to find the optimal pairwise (k and $\lambda$)
::
  % Snakefile
  % snakemake -j num_cpu 

:: Step 2: In the Sanakefile, provide the range of the number of signatures $k$ and $lambda$ for the grid search to determine the optimal k and lambda.  

..code-block::bash
  localrules: all
     ks = list(range(3, 6))
     lambdas = [0, 0.01,0.05, 0.1, 0.2]


Contact
-------

Maintainer Name - chrisbop@uio.no
