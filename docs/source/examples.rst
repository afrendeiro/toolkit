Examples
******************************


Analysis example
==============================

The following is an example of how to use ``ngs_toolkit`` in a ATAC-seq project.
While straightforward, it still allows considerable customization due to the modularity of the toolkit and the parametrization of most functions (this example uses default values everywhere nonetheless).


We have the following `PEP project <https://peppy.readthedocs.io>`_ config YAML file:

.. code-block:: yaml

    project_name: example_project
    project_description: example_project
    username: user
    email: user@cemm.oeaw.ac.at
    metadata:
      output_dir: /scratch/lab_bock/shared/projects/example_project
      results_subdir: data
      submission_subdir: submission
      pipeline_interfaces: /home/user/workspace/open_pipelines/pipeline_interface.yaml
      sample_annotation: /scratch/lab_bock/shared/projects/example_project/metadata/annotation.csv
      sample_subannotation: /scratch/lab_bock/shared/projects/example_project/metadata/sample_subannotation.csv
      comparison_table: /scratch/lab_bock/shared/projects/example_project/metadata/comparison_table.csv
    data_sources:
      bsf: /path/to/samples/{flowcell}/{flowcell}_{lane}#{sample_name}.bam
    genomes:
      human: hg19
    trackhubs:
      trackhub_dir: /path/to/public_html/user/example_project/
      url: http://root-url.com/example_project


The following sample annotation CSV file:

.. csv-table:: Annotation table for example
   :header: "sample_name", "genotype", "replicate", "organism", flowcell, lane

    "ATAC-seq_KOA_r1",  "KO_A",   "1",   "human", "C0RQ31ACXXX",   "1"
    "ATAC-seq_KOA_r2",  "KO_A",   "2",   "human", "C0RQ31ACXXX",   "1"
    "ATAC-seq_KOB_r1",  "KO_B",   "1",   "human", "C0RQ31ACXXX",   "1"
    "ATAC-seq_KOB_r2",  "KO_B",   "2",   "human", "C0RQ31ACXXX",   "1"
    "ATAC-seq_WT_r1",   "WT",   "1",    "human",    "C0RQ31ACXXX", "1"
    "ATAC-seq_WT_r2",   "WT",    "2",   "human", "C0RQ31ACXXX",    "1"


And the following comparison table:

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
    analysis = ATACSeqAnalysis(
        from_pep=os.path.join("metadata", "project_config.yaml"))

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

    # Save analysis object
    analysis.to_pickle()


    # UNSUPERVISED ANALYSIS
    # # plot pairwise sample correlations, 
    # # perform dimensionality reduction (MDS, PCA)
    # # and plot samples in this spaces, annotated with their attributes
    analysis.unsupervised_analysis()


    # SUPERVISED ANALYSIS
    # # differential analysis with DESeq2
    analysis.differential_analysis()

    # # Save analysis object
    analysis.to_pickle()

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
