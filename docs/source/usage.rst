Quick usage
=============================


Interactive usage through the API
-----------------------------

To use a particular class or function from the toolkit, simply import it
following the structure of the library:

.. code-block:: python

    from ngs_toolkit import ATACSeqAnalysis
    from ngs_toolkit.utils import log_p_values

The :class:`ngs_toolkit.analysis.Analysis` and their data type-specific
children are the main drivers of the workflow, storing attributes and providing
various methods through an OOP interface:

.. code-block:: python

    from ngs_toolkit.demo import generate_project

    an = generate_project(data_type="ATAC-seq", sample_input_files=True)
    an.measure_coverage()
    an.normalize()
    an.unsupervised_analysis()
    an.differential_analysis()
    an.plot_differential()
    an.get_peak_gene_annotation()
    an.annotate_features()
    an.differential_enrichment(steps=['enrichr'])
    an.plot_differential_enrichment()


Running recipes through the command-line interface
-----------------------------

``ngs_toolkit`` also has some command-line programs on some commonly used
workflows (here called ``recipes``), which can be run in the following manner:

.. code-block:: bash

    PEP=`python -m ngs_toolkit.recipes.generate_project --sample-input-files True`
    python -m ngs_toolkit.recipes.ngs_analysis $PEP

This example is roughly equivalent to the on above with interactive usage.
