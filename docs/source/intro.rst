Introduction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Installation
=============================

``ngs_toolkit`` is available for Python 3 only.

To install, simply do:

.. code-block:: bash

   pip install ngs-toolkit [--user]

you might need the ``--user`` flag if not root or running in a virtual environment.

This will install all the Python dependencies needed too. See `here <https://github.com/afrendeiro/toolkit/blob/master/requirements/requirements.txt>`_ a list of all Python dependencies used.


If you wish to install libraries required for additional work, you can do:

.. code-block:: bash

   pip install ngs-toolkit[rstats] [--user]


Non-Python optional requirements
-----------------------------

``ngs_toolkit`` makes use of some non-Python dependencies.

The following are either required only for some data types or just recommended:

 - `bedtools <https://bedtools.readthedocs.io/en/latest/>`_: required for some ATAC/ChIP-seq functions. It is underlying the Python interface library to bedtools (pybedtools) which can be installed without bedtools.
 - `R <https://www.r-project.org/>`_ and some bioconductor libraries (optional):

   - `cqn <https://bioconductor.org/packages/release/bioc/html/cqn.html>`_: used for GC-content aware normalization of NGS data.
   - `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_: used for differential testing of genes/regulatory elements.
 - `Kent tools <https://github.com/ENCODE-DCC/kentUtils>`_ (optional): the '2bitToFa' binary from UCSC's Kent bioinformatics toolkit is used to convert between the 2bit and FASTA formats.


API usage
=============================

To use a particular class or function from the toolkit, import it like this from within Python/iPython:

.. code-block:: python

   from ngs_toolkit import ATACSeqAnalysis
   from ngs_toolkit.utils import log_pvalues
