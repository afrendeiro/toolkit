Introduction
^^^^^^^^

Installation
------------

To install ``ngs_toolkit`` simply do:

.. code-block:: bash

   pip install https://github.com/afrendeiro/toolkit/zipball/master [--user]

you might need the ``--user`` flag if not root or running in a virtual environment.


API usage
---------

To use a particular class or function from the toolkit, import it like this:

.. code-block:: python

   from ngs_toolkit.atacseq import ATACSeqAnalysis
   

Full analysis example usage
---------------------------

The following is a real-world example of how to use ``ngs_toolkit`` in a ATAC-seq project.
While straightforward, it still allows considerable customization due to the modularity of the toolkit and the parametrization of most functions (this example uses default values everywhere nonetheless).

.. code-block:: python

    import os
    from looper.models import Project
    from ngs_toolkit.atacseq import ATACSeqAnalysis

    # Start project and analysis objects
    prj = Project(os.path.join("metadata", "project_config.yaml"))
    prj.samples = [sample for sample in prj.samples if sample.library == "ATAC-seq"]
    atac_analysis = ATACSeqAnalysis(name="example_analysis", prj=prj, samples=prj.samples)

    ## Project's attributes to use in analysis
    sample_attributes = ["sample_name", "library", "patient_id", "visit", "cell_number", "sex", "batch", "flowcell", "lane"]
    attributes_to_plot = ["patient_id", "visit", "cell_number", "sex", "batch", "flowcell", "lane"]

    # Generate consensus peak set and annotate it
    ## get consensus peak set from all samples
    atac_analysis.get_consensus_sites()
    ## calculate peak support (what fraction of samples has which peak)
    atac_analysis.calculate_peak_support()
    ## annotate peak set with genes
    atac_analysis.get_peak_gene_annotation()
    ## annotate peak set with genomic context
    atac_analysis.get_peak_genomic_location()
    ## annotate peak set with chromatin context
    atac_analysis.get_peak_chromatin_state(
        os.path.join(atac_analysis.data_dir, "external", "E032_15_coreMarks_mnemonics.bed"))

    # Use accessibility quantitatively
    ## get coverage values for each peak in each sample of ATAC-seq
    atac_analysis.measure_coverage()

    # Normalize accessibility (quantile normalization + GC correction)
    atac_analysis.normalize(method="gc_content")

    # Annotate normalized accessibility with sample and region info
    # annotate matrix with peak metadata
    atac_analysis.annotate()
    # annotate matrix with sample metadata
    atac_analysis.accessibility = atac_analysis.annotate_with_sample_metadata(
        quant_matrix="coverage_annotated",
        attributes=sample_attributes)

    # Save analysis object
    atac_analysis.to_pickle()

    # Perform unsupervised analysis
    atac_analysis.unsupervised(
        quant_matrix="accessibility", samples=None,
        attributes_to_plot=attributes_to_plot, plot_prefix="accessibility")


    # The analysis can continue with more supervised methods (appropriate to the project at hand) with the following methods:

    from ngs_toolkit.general import (differential_analysis,
                                 plot_differential,
                                 differential_enrichment,
                                 collect_differential_enrichment,
                                 plot_differential_enrichment)
