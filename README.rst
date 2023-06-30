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
  
** Running mutational signature refit algorithm:**
::
  $cucumber refit --help
:: Running cucumber refit using mutational catalogue matrix
  $cucumber refit example_data/skin20.txt example_data/sig_cosmic_v3_2019.txt output/expo.txt output/expo_avg.txt


Contact
-------

Maintainer Name - chrisbop@uio.no
