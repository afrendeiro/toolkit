Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With pip
=============================

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
=============================

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


Non-Python requirements
=============================


``ngs_toolkit`` makes use of some non-Python dependencies.

 - `bedtools <https://bedtools.readthedocs.io/en/latest/>`_: version should be at least 2.27.1

The following are highly recommended only for some data or analysis types:

 - `R <https://www.r-project.org/>`_ and some bioconductor libraries (optional):
   - `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_ (optional): used for differential testing of genes/regulatory elements and variance stabilization transformation of data.
   - `cqn <https://bioconductor.org/packages/release/bioc/html/cqn.html>`_ (optional): used for GC-content aware normalization of NGS data.
 - `Kent tools <https://github.com/ENCODE-DCC/kentUtils>`_ (optional):
   - the ``twoBitToFa`` binary from UCSC's Kent bioinformatics toolkit is used to convert between the 2bit and FASTA formats.

For region-based enrichment analysis, you may also want to have the following software installed (entirely optional):

 - `MEME suite <http://meme-suite.org/>`_
 - `HOMER motif analysis <http://homer.ucsd.edu/homer/motif/>`_
 - `LOLA R package <http://code.databio.org/LOLA/>`_

You can see how to install all requirements in an Ubuntu-based system in the provided `Dockerfile <https://github.com/afrendeiro/toolkit/blob/master/Dockerfile>`_.


Docker
=============================

A Docker image containing ``ngs_toolkit`` and its dependencies is also available: https://hub.docker.com/r/afrendeiro/ngs-toolkit

To pull the image and run a module for example in this way:

.. code-block:: bash

   docker pull afrendeiro/ngs-toolkit
   docker run ngs-toolkit python3 -m ngs_toolkit.recipes.ngs_analysis --help

You can also run an interactive session of ``ngs_toolkit`` `based on the docker image on Gitpod <https://gitpod.io/#https://github.com/afrendeiro/toolkit>`_.

The Dockerfile that produced the image is available in the github repository: https://github.com/afrendeiro/toolkit/blob/master/Dockerfile
