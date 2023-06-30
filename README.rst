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
   $ pip install cucumber
To obtain help :
::
   $ cucumber --help
   $ Usage: cucumber [OPTIONS] COMMAND [ARGS]...
Commands:
::
  count-mutation  Count mutation types in a VCF file
  denovo Parameters ---matrix_file: str n_signatures 
  refit  Parameters ---numeric_chromosomes n_bootstraps
  
Running mutational signature refit algorithm:
::
$cucumber refit --help
Usage: cucumber refit [OPTIONS] MATRIX_FILE SIGNATURE_FILE
                      OUTPUT_FILE_EXPOSURE OUTPUT_FILE_EXPOSURE_AVG

  Parameters --- numeric_chromosomes n_bootstraps matrix_file: str
  signature_file: str output_file_exposure: str opportunity_file: str
  data_type: DataType numeric_chromosomes: bool     True if chromosome names
  in vcf are '1', '2', '3'. False if 'chr1', 'chr2', 'chr3' genotyped: bool
  True if the VCF file has genotype information for many samples

Arguments:
  MATRIX_FILE               [required]
  SIGNATURE_FILE            [required]
  OUTPUT_FILE_EXPOSURE      [required]
  OUTPUT_FILE_EXPOSURE_AVG  [required]

Options:
  --opportunity-file TEXT
  --data-type [exome|genome]      [default: DataType.exome]
  --n-bootstraps INTEGER          [default: 100]
  --numeric-chromosomes / --no-numeric-chromosomes
                                  [default: no-numeric-chromosomes]
  --genotyped / --no-genotyped    [default: genotyped]
  --output-file-exposure-plot TEXT
  --help                          Show this message and exit.

Contact
-------

Maintainer Name - chrisbop@uio.no
