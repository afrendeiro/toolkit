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
      sample_subannotation: /scratch/lab_bock/shared/projects/example_project/metadata/merge_table.csv
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



ATAC-seq example
-------------------------------

.. code-block:: python

    import os
    from looper.models import Project
    from ngs_toolkit.atacseq import ATACSeqAnalysis    
    from ngs_toolkit.general import (differential_analysis,
                                 plot_differential,
                                 differential_enrichment,
                                 plot_differential_enrichment)


    # Start project and analysis objects
    prj = Project(
        os.path.join("metadata", "project_config.yaml"))
    prj.samples = [sample for sample in prj.samples if sample.library == "ATAC-seq"]
    atac_analysis = ATACSeqAnalysis(
        name="example.atac_analysis",
        prj=prj,
        samples=prj.samples)
    data_type = "ATAC-seq"

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
        os.path.join(
            atac_analysis.data_dir,
            "external",
            "E032_15_coreMarks_mnemonics.bed"))

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
        attributes=prj.sample_attributes)

    # Save analysis object
    atac_analysis.to_pickle()


    # UNSUPERVISED ANALYSIS

    # plot pairwise sample correlations, 
    # perform dimensionality reduction (MDS, PCA)
    # and plot samples in this spaces, annotated with their attributes
    atac_analysis.unsupervised(
        quant_matrix="accessibility", samples=None,
        attributes_to_plot=attributes_to_plot, plot_prefix="accessibility")


    # SUPERVISED ANALYSIS

    # read in comparison table, subset if needed
    comparison_table = pd.read_csv(os.path.join("metadata", "comparison_table.csv"))
    comparison_table = comparison_table[
        (comparison_table['data_type'] == data_type) &
        (comparison_table['comparison_type'] == 'differential')]

    # differential analysis with DESeq2
    analysis.differential_results = differential_analysis(
        analysis,
        comparison_table,
        data_type=data_type,
        output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
        covariates=None,
        overwrite=True)

    # Save analysis object
    analysis.to_pickle()

    # plot scatter, volcano, MA, heatmaps on the differential regions
    # by groups and with individual samples, with normalized values
    # and scalled values (Z-score).
    plot_differential(
        analysis,
        analysis.differential_results,
        matrix=getattr(analysis, quant_matrix),
        comparison_table=comparison_table,
        output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
        output_prefix="differential_analysis",
        data_type=data_type,
        alpha=alpha,
        corrected_p_value=True,
        fold_change=abs_fold_change,
        rasterized=True,
        robust=True,
        group_wise_colours=True,
        group_variables=group_variables)

    # perform enrichment analysis on differnetial region sets
    # using LOLA, MEME-AME, HOMER and Enrichr
    differential_enrichment(
        analysis,
        diff,
        data_type=data_type,
        output_dir=output_dir,
        directional=True,
        max_diff=n_top,
        sort_var="pvalue",
        as_jobs=False)

    # for each type of enrichment results,
    # plot bar and scatter plots of odds ratio vs p-value,
    # heatmaps of enrichment across terms for each comparison
    # and comparison correlation in enrichment terms
    for enrichment_name, enrichment_type in [
            ('motif', 'meme_ame'), ('homer_consensus', 'homer_consensus'),
            ('lola', 'lola'), ('enrichr', 'enrichr')]:
        enrichment_table = pd.read_csv(
                os.path.join(output_dir, "differential_analysis" + ".{}.csv".format(enrichment_type)))

        plot_differential_enrichment(
            enrichment_table,
            enrichment_name,
            data_type=data_type,
            output_dir=output_dir,
            direction_dependent=True, barplots=False,
            top_n=10 if enrichment_name not in ["motif", "homer_consensus"] else 50)
