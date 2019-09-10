Introduction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Installation
=============================

With pip
-----------------------------

``ngs_toolkit`` is available for Python 3 only. Is is tested in Python 3.6 and 3.7.

To install, simply do:

.. code-block:: bash

   pip install ngs-toolkit

you might need to add a ``--user`` flag if not root or running in a virtual environment.

This will install all the Python dependencies needed too.
See `here <https://github.com/afrendeiro/toolkit/blob/master/requirements/requirements.txt>`_ a list of all Python dependencies used.

To install the latest development version:

.. code-block:: bash

   pip install git+https://github.com/afrendeiro/toolkit.git#egg=ngs-toolkit


Using a conda environment
-----------------------------

Get the `latest Python 3 installation of miniconda from the conda website <https://docs.conda.io/en/latest/miniconda.html>`_ and follow the instructions for installation and activation of the environment.

Setup the bioconda channel:

.. code-block:: bash

	conda config --add channels defaults
	conda config --add channels bioconda
	conda config --add channels conda-forge

Install non-Python dependencies:

.. code-block:: bash

   conda install -y bedtools==2.27.1
   conda install -y ucsc-twobittofa
   conda install -y bioconductor-deseq2
   conda install -y bioconductor-cqn

And then install the ``ngs-toolkit`` library with pip (available only through PyPi).

.. code-block:: bash

	pip install ngs-toolkit


**Non-Python requirements**
-----------------------------


``ngs_toolkit`` makes use of some non-Python dependencies.

 - `bedtools <https://bedtools.readthedocs.io/en/latest/>`_: version should be at least 2.27.1

The following are highly recommended only for some data or analysis types:

 - `R <https://www.r-project.org/>`_ and some bioconductor libraries (optional):
   - `cqn <https://bioconductor.org/packages/release/bioc/html/cqn.html>`_ (optional): used for GC-content aware normalization of NGS data.
   - `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_ (optional): used for differential testing of genes/regulatory elements.
 - `Kent tools <https://github.com/ENCODE-DCC/kentUtils>`_ (optional): the 'twoBitToFa' binary from UCSC's Kent bioinformatics toolkit is used to convert between the 2bit and FASTA formats.


Testing
=============================

To make sure everything is correctly configured, the user is encouraged to test the library prior to use.

In order to do this, install testing requirements and simply run ``pytest``:

.. code-block:: bash

   pip install ngs-toolkit[testing]
   pytest --pyargs ngs_toolkit


Pytest will output summary results (`see for example <https://travis-ci.org/afrendeiro/toolkit/jobs/580167563>`_) and further outputs can be seen in ``${TMPDIR}/pytest-of-${USER}/`` or ``/tmp/pytest-of-${USER}/`` if $TMPDIR is not defined.


API usage
=============================

To use a particular class or function from the toolkit, import it like this from within Python/iPython:

.. code-block:: python

   from ngs_toolkit import ATACSeqAnalysis
   from ngs_toolkit.utils import log_pvalues
