Examples
******************************


Analysis example
==============================

The following is an example of how to use ``ngs_toolkit`` in a ATAC-seq project.
While straightforward, it still allows considerable customization due to the modularity of the toolkit and the parametrization of most functions (this example uses default values everywhere nonetheless).

.. note::
    ``ngs_toolkit`` from version 0.25.0 on uses the `PEP 2.0 specification <http://pep.databio.org/en/2.0.0/specification/>`_. If you have a PEP made in an earlier version, you must update it in order to use ``ngs_toolkit``>=0.25.0.


We have the following `PEP project <http://pep.databio.org/>`_ config YAML file:

.. code-block:: yaml

    pep_version: "2.0.0"
    name: example_project
    description: example_project
    username: user
    email: user@email.com

    sample_table: annotation.csv
    subsample_table:
    comparison_table: comparison_table.csv

    submission_subdir: submission
    results_subdir: data
    output_dir: example_project

    pipeline_interfaces: /home/user/workspace/open_pipelines/pipeline_interface.yaml

    sample_attributes:
      - sample_name
      - genotype
      - replicate
    group_attributes:
      - genotype
      - replicate
    sample_modifiers:
        imply:
            - if:
                organism: 'human'
              then:
                genome: 'hg38'
        derive:
            attributes: [data_source]
            sources:
                local: data/{sample_name}.bam
                bsf: /scratch/lab_bsf/samples/{flowcell}/{flowcell}_{lane}_samples/{flowcell}_{lane}#{BSF_name}.bam



The following sample annotation CSV file, 'annotation.csv':

.. csv-table:: Annotation table for example
   :header: "sample_name", "protocol", "genotype", "replicate", "organism", flowcell, lane

    "ATAC-seq_KOA_r1",  "ATAC-seq", KO_A",   "1",   "human", "C0AXX",   "1"
    "ATAC-seq_KOA_r2",  "ATAC-seq", KO_A",   "2",   "human", "C0AXX",   "1"
    "ATAC-seq_KOB_r1",  "ATAC-seq", KO_B",   "1",   "human", "C0AXX",   "1"
    "ATAC-seq_KOB_r2",  "ATAC-seq", KO_B",   "2",   "human", "C0AXX",   "1"
    "ATAC-seq_WT_r1",   "ATAC-seq", WT",     "1",   "human", "C0AXX",   "1"
    "ATAC-seq_WT_r2",   "ATAC-seq", WT",     "2",   "human", "C0AXX",   "1"


And the following comparison table, 'comparison_table.csv':

.. csv-table:: Comparison table for example
   :header: "comparison_name", "comparison_side", "sample_name", "sample_group"
   :widths: 30, 30, 30, 30

    "KOA_vs_WT",    "1",    "ATAC-seq_KOA_r1",  "KO_A"
    "KOA_vs_WT",    "1",    "ATAC-seq_KOA_r2",  "KO_A"
    "KOA_vs_WT",    "0",    "ATAC-seq_WT_r1",   "WT"
    "KOA_vs_WT",    "0",    "ATAC-seq_WT_r2",   "WT"
    "KOB_vs_WT",    "1",    "ATAC-seq_KOB_r1",  "KO_B"
    "KOB_vs_WT",    "1",    "ATAC-seq_KOB_r2",  "KO_B"
    "KOB_vs_WT",    "0",    "ATAC-seq_WT_r1",   "WT"
    "KOB_vs_WT",    "0",    "ATAC-seq_WT_r2",   "WT"



ATAC-seq analysis example
-------------------------------

.. code-block:: python

    import os
    from ngs_toolkit.atacseq import ATACSeqAnalysis

    # Start project and analysis objects
    analysis = ATACSeqAnalysis(from_pep="project_config.yaml")

    # Generate consensus peak set and annotate it
    ## get consensus peak set from all samples
    analysis.get_consensus_sites()
    ## annotate peak set with genomic context
    analysis.get_peak_genomic_location()
    ## annotate peak set with chromatin context
    analysis.get_peak_chromatin_state(
        os.path.join(
            analysis.data_dir,
            "external",
            "E032_15_coreMarks_mnemonics.bed"))
    ## annotate peak set with genes
    analysis.get_peak_gene_annotation()

    # Use accessibility quantitatively
    ## get coverage values for each peak in each sample of ATAC-seq
    analysis.measure_coverage()

    # Normalize accessibility (quantile normalization + GC correction, requires cqn R library)
    analysis.normalize(method="cqn")

    # Annotate normalized accessibility with sample and region info
    # # annotate dataframe with peak metadata
    analysis.annotate_features()
    # # annotate dataframe with sample metadata
    analysis.accessibility = analysis.annotate_samples()

    # UNSUPERVISED ANALYSIS
    # # plot pairwise sample correlations,
    # # perform dimensionality reduction (MDS, PCA)
    # # and plot samples in this spaces, annotated with their attributes
    analysis.unsupervised_analysis()


    # SUPERVISED ANALYSIS
    # # differential analysis with DESeq2
    analysis.differential_analysis()

    # # plot scatter, volcano, MA, heatmaps on the differential regions
    # # by groups and with individual samples, with normalized values
    # # and scalled values (Z-score).
    analysis.plot_differential(
        alpha=0.05,
        corrected_p_value=True,
        fold_change=1)

    # # perform enrichment analysis on differnetial region sets
    # # using LOLA, MEME-AME, HOMER and Enrichr
    analysis.differential_enrichment(
        directional=True,
        max_diff=1000,
        sort_var="pvalue")

    # # for each type of enrichment results,
    # # plot bar and scatter plots of odds ratio vs p-value,
    # # heatmaps of enrichment across terms for each comparison
    # # and comparison correlation in enrichment terms
    analysis.plot_differential_enrichment()
