Introduction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Installation
=============================

With pip
-----------------------------

``ngs_toolkit`` is available for Python 3 only.

To install, simply do:

.. code-block:: bash

   pip install ngs-toolkit

you might need to add a ``--user`` flag if not root or running in a virtual environment.

This will install all the Python dependencies needed too.
See `here <https://github.com/afrendeiro/toolkit/blob/master/requirements/requirements.txt>`_ a list of all Python dependencies used.

If you wish to install optional libraries that interface with R libraries, you can pass ``[rstats]`` to the following pip call:

.. code-block:: bash

   pip install ngs-toolkit[rstats]


To install the latest development version:

.. code-block:: bash

   pip install git+https://github.com/afrendeiro/toolkit.git#egg=ngs-toolkit



**Non-Python optional requirements**


``ngs_toolkit`` makes use of some non-Python dependencies.

The following are required only for some data or analysis types:

 - `bedtools <https://bedtools.readthedocs.io/en/latest/>`_: required for some ATAC/ChIP-seq functions. It is underlying the Python interface library to bedtools (pybedtools) which can be installed without bedtools.
 - `R <https://www.r-project.org/>`_ and some bioconductor libraries (optional):

   - `cqn <https://bioconductor.org/packages/release/bioc/html/cqn.html>`_: used for GC-content aware normalization of NGS data.
   - `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_: used for differential testing of genes/regulatory elements.
 - `Kent tools <https://github.com/ENCODE-DCC/kentUtils>`_ (optional): the 'twoBitToFa' binary from UCSC's Kent bioinformatics toolkit is used to convert between the 2bit and FASTA formats.

.. note::
	``bedtools`` version should be below 2.24.0 (2.20.1 is used for testing)


Using a conda environment
-----------------------------

Get the `latest Python 3 installation of miniconda from the conda website <https://docs.conda.io/en/latest/miniconda.html>`_ and follow the instructions for installation and activation of the environment.

Setup the bioconda channel:

.. code-block:: bash

	conda config --add channels defaults
	conda config --add channels bioconda
	conda config --add channels conda-forge

Install the dependencies:

.. code-block:: bash

   conda install -y bedtools==2.20.1
   conda install -y ucsc-twobittofa
   conda install -y bioconductor-deseq2
   conda install -y bioconductor-cqn

And then install the ``ngs-toolkit`` library with pip (available only through PyPi).

.. code-block:: bash

	pip install ngs-toolkit


API usage
=============================

To use a particular class or function from the toolkit, import it like this from within Python/iPython:

.. code-block:: python

   from ngs_toolkit import ATACSeqAnalysis
   from ngs_toolkit.utils import log_pvalues
