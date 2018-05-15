Examples
******************************


ATAC-seq analysis example
---------------------------

The following is an example of how to use ``ngs_toolkit`` in a ATAC-seq project.
While straightforward, it still allows considerable customization due to the modularity of the toolkit and the parametrization of most functions (this example uses default values everywhere nonetheless).


We have the following PEP/looper project config YAML file:

.. code-block:: yaml

    project_name: baf_complex
    project_description: baf_complex
    username: user
    email: user@cemm.oeaw.ac.at
    metadata:
      output_dir: /scratch/lab_bock/shared/projects/baf_complex
      results_subdir: data
      submission_subdir: submission
      pipeline_interfaces: /home/user/workspace/open_pipelines/pipeline_interface.yaml
      sample_annotation: /scratch/lab_bock/shared/projects/baf_complex/metadata/annotation.csv
      merge_table: /scratch/lab_bock/shared/projects/baf_complex/metadata/merge_table.csv
    data_sources:
      bsf: /path/to/samples/{flowcell}/{flowcell}_{lane}#{sample_name}.bam
    genomes:
      human: hg19
    trackhubs:
      trackhub_dir: /data/groups/lab_bock/public_html/user/baf_complex/
      url: http://biomedical-sequencing.at/bocklab/user/baf_complex


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


Example analysis:


.. code-block:: python

    import os
    from looper.models import Project
    from ngs_toolkit.atacseq import ATACSeqAnalysis    
    from ngs_toolkit.general import (differential_analysis,
                                 plot_differential,
                                 differential_enrichment,
                                 collect_differential_enrichment,
                                 plot_differential_enrichment)


    # Start project and analysis objects
    prj = Project(
        os.path.join("metadata", "project_config.yaml"))
    prj.samples = [sample for sample in prj.samples if sample.library == "ATAC-seq"]
    atac_analysis = ATACSeqAnalysis(
        name="example.atac_analysis",
        prj=prj,
        samples=prj.samples)

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
    # (all at once)
    try:
        analysis.differential_results = differential_analysis(
            analysis,
            comparison_table,
            data_type=data_type,
            samples=[s for s in analysis.samples if s.name in comparison_table['sample_name'].tolist()],
            output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
            covariates=None,
            alpha=alpha,  # not really used 
            overwrite=True)
    # (one comparison at a time)
    except:
        analysis.differential_results = pd.DataFrame()
        for comparison in comparison_table['comparison_name'].unique():
            comp = comparison_table[comparison_table['comparison_name'] == comparison_name]
            res = differential_analysis(
                analysis,
                comp,
                data_type=data_type,
                samples=[s for s in analysis.samples if s.name in comp['sample_name'].tolist()],
                output_dir="{}/differential_analysis_{}".format(analysis.results_dir, data_type),
                covariates=None,
                alpha=alpha,  # not really used 
                overwrite=True)
            analysis.differential_results = analysis.differential_results.append(res, ignore_index=True)

    analysis.differential_results = analysis.differential_results.set_index("index")
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
