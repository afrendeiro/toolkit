#!/usr/bin/env python

"""
This is "ngs_analysis" recipe for ngs_toolkit.
"""


if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')


from argparse import ArgumentParser
import os
import sys

import matplotlib
import pandas as pd
import seaborn as sns

from peppy import Project
from ngs_toolkit.atacseq import ATACSeqAnalysis
from ngs_toolkit.rnaseq import RNASeqAnalysis
from ngs_toolkit.chipseq import ChIPSeqAnalysis
from ngs_toolkit.general import (collect_differential_enrichment,
                                 differential_analysis,
                                 differential_enrichment, differential_overlap,
                                 plot_differential,
                                 plot_differential_enrichment,
                                 unsupervised_analysis)


# Set settings
pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


class Unbuffered(object):
    def __init__(self, stream):
        self.stream = stream

    def write(self, data):
        self.stream.write(data)
        self.stream.flush()

    def writelines(self, datas):
        self.stream.writelines(datas)
        self.stream.flush()

    def __getattr__(self, attr):
        return getattr(self.stream, attr)


sys.stdout = Unbuffered(sys.stdout)


def add_args(parser):
    """
    Global options for analysis.
    """
    parser.add_argument(
        dest="config_file",
        help="YAML project configuration file.",
        type=str)
    parser.add_argument(
        "-n", "--analysis-name",
        dest="name",
        default=None,
        help="Name of analysis. Will be the prefix of output_files. "
             "By default it will be the name of the Project given in the YAML configuration.",
        type=str)
    parser.add_argument(
        "-o", "--results-output",
        default="results",
        dest="results_dir",
        help="Directory for analysis output files. "
        "Default is 'results' under the project roort directory.",
        type=str)
    parser.add_argument(
        "-t", "--data-type",
        default=None,
        choices=["ATAC-seq", "RNA-seq", "ChIP-seq"],
        dest="data_type",
        help="Data type to restrict analysis to. "
        "Default is to run separate analysis for each data type.",
        type=str)
    parser.add_argument(
        "-q", "--pass-qc",
        action="store_true",
        dest="pass_qc",
        help="Whether only samples with a 'pass_qc' value of '1' "
        "in the annotation sheet should be used.")
    parser.add_argument(
        "-a", "--alpha",
        default=0.05,
        dest="alpha",
        help="Alpha value of confidence for supervised analysis.",
        type=str)
    parser.add_argument(
        "-f", "--fold-change",
        default=0,
        dest="abs_fold_change",
        help="Absolute log2 fold change value for supervised analysis.",
        type=str)
    return parser


def main():
    parser = ArgumentParser(
        prog="ngs_analysis_recipe",
        description="NGS analysis recipe."
    )
    parser = add_args(parser)
    args = parser.parse_args()
    # args = parser.parse_args('-t ATAC-seq metadata/project_config.yaml'.split(" "))

    # Start project
    print("Starting peppy project with project configuration file: '{}'".format(args.config_file))
    prj = Project(args.config_file)
    print("Changing directory to project root directory: '{}'.".format(prj.metadata.output_dir))
    os.chdir(prj.metadata.output_dir)
    if args.pass_qc:
        print("Filtering samples out which didn't pass QC as specified in sample annotation in column 'pass_qc'")
        prj._samples = [s for s in prj._samples if s.pass_qc not in ['0', 0, 'False', False]]
    print("Setting location of sample files dependent on sample types.")
    for sample in prj.samples:
        if hasattr(sample, "protocol"):
            sample.library = sample.protocol

        if sample.library in ["ATAC-seq", "ChIP-seq", "ChIPmentation"]:
            sample.mapped = os.path.join(
                sample.paths.sample_root,
                "mapped", sample.name + ".trimmed.bowtie2.bam")
            sample.filtered = os.path.join(
                sample.paths.sample_root,
                "mapped", sample.name + ".trimmed.bowtie2.filtered.bam")
            sample.peaks = os.path.join(
                sample.paths.sample_root,
                "peaks", sample.name + "_peaks.narrowPeak")
        elif sample.library == "RNA-seq":
            sample.bitseq_counts = os.path.join(
                sample.paths.sample_root, "bowtie1_{}".format(sample.transcriptome),
                "bitSeq", sample.name + ".counts")

    # ANALYSIS
    if args.data_type is None:
        print("Type of analysis not specified. Will run independent analysis for all types of data in the sample annotation sheet.")
        data_types = sorted(list(set([s.library for s in prj.samples])))
        print("Sample data types: '{}'.".format(",".join(data_types)))
    else:
        print("Type of analysis specified. Will run only analysis for samples of type '{}'.".format(args.data_type))
        data_types = [args.data_type]
        print("Sample data types: '{}'.".format(",".join(data_types)))
    if args.name is None:
        print("Analysis name not specified, will use name in project configuration file: '{}'.".format(prj.project_name))
        args.name = prj.project_name

    for data_type in data_types:
        print("Starting analysis for samples of type: '{}'.".format(data_type))
        samples = [s for s in prj.samples if (s.library == data_type)]
        if len(samples) > 0:
            print(
                "Samples under consideration: '{}'. ".format(",".join([s.name for s in samples])) +
                "Total of {} samples.".format(len([s.name for s in samples])))
        else:
            raise ValueError("There were no valid samples for this analysis type!")

        if data_type in ["ATAC-seq"]:
            print("Initializing ATAC-seq analysis")
            analysis = ATACSeqAnalysis(
                name=args.name + "_atacseq", prj=prj,
                samples=samples, results_dir=args.results_dir)
        elif data_type in ["ChIP-seq"]:
            print("Initializing ChIP-seq analysis")
            analysis = ChIPSeqAnalysis(
                name=args.name + "_chipseq", prj=prj,
                samples=samples, results_dir=args.results_dir)
        elif data_type in ["RNA-seq"]:
            print("Initializing RNA-seq analysis")
            analysis = RNASeqAnalysis(
                name=args.name + "_rnaseq", prj=prj,
                samples=samples, results_dir=args.results_dir)

        if hasattr(prj, "sample_attributes"):
            print("Using sample attributes from project configuration file: '{}'.".format(",".join(prj.sample_attributes)))
            sample_attributes = prj.sample_attributes
        else:
            print("Project configuration file does not contain a 'sample_attributes' section.")
            print("Sample annotation will be minimal!")
            sample_attributes = ['sample_name']
        if hasattr(prj, "group_attributes"):
            print("Using group attributes from project configuration file: '{}'.".format(",".join(prj.group_attributes)))
            group_attributes = prj.group_attributes
        else:
            print("Project configuration file does not contain a 'group_attributes' section.")
            print("Group-wise labeling of samples will not be possible!")
            group_attributes = ['sample_name']

        if "comparison_table" in prj.metadata.keys():
            comparison_table_file = prj.metadata['comparison_table']
            print("Using comparison table specified in project configuration file: '{}'.".format(comparison_table_file))
        else:
            comparison_table_file = os.path.join("metadata", "comparison_table.csv")
            print("Will try using comparison table in metadata directory: '{}'.".format(comparison_table_file))

        try:
            comparison_table = pd.read_csv(comparison_table_file)
        except IOError:
            print("Comparison table not present for project.")
            print("Will not perform comparisons.")
            comparison_table = pd.DataFrame()

        return main_analysis_pipeline(
            analysis, data_type=data_type,
            sample_attributes=sample_attributes,
            plotting_attributes=group_attributes,
            comparison_table=comparison_table,
            alpha=args.alpha, abs_fold_change=args.abs_fold_change
        )


def main_analysis_pipeline(
        analysis, data_type,
        sample_attributes, plotting_attributes,
        comparison_table,
        alpha=0.05, abs_fold_change=0):
    """
    Main analysis pipeline for ATAC-seq and RNA-seq - untested for ChIP-seq/ChIPmentation.
    Gets quantification matrices, normalizes them,
    performes unsupervised and supervised analysis and
    gets and plots enrichments for supervised analysis.
    """

    genomes = list(set(s.genome for s in analysis.samples))

    if len(genomes) != 1:
        raise ValueError("Samples under analysis have more than one genome assembly: '{}'.".format("', '".join(genomes)))
    else:
        genome = genomes[0]

    if data_type in ["ATAC-seq", "ChIP-seq"]:

        # GET CONSENSUS PEAK SET, ANNOTATE IT, PLOT
        # Get consensus peak set from all samples
        if data_type == "ChIP-seq":
            comps = comparison_table[
                (comparison_table['data_type'] == data_type) &
                (comparison_table['comparison_type'] == 'peaks')]
            if comps.empty:
                print("Comparison table has no comparisons with 'data_type'=='ChIP-seq' and 'comparison_type'=='peaks'.")
                print("Aborting.")
                return 1
            analysis.get_consensus_sites(comps, region_type="summits")
        else:
            analysis.get_consensus_sites(region_type="summits")
        analysis.calculate_peak_support(region_type="summits")

        # GET CHROMATIN OPENNESS MEASUREMENTS, PLOT
        # Get coverage values for each peak in each sample of ATAC-seq and ChIPmentation
        analysis.measure_coverage()
        # normalize coverage values
        analysis.normalize_coverage_quantiles()
        analysis.get_peak_gccontent_length(genome=genome)
        analysis.normalize_gc_content()

        # Annotate peaks with closest gene
        if os.path.exists(os.path.join("data", "external", "refseq.refflat.tss.bed")):
            analysis.get_peak_gene_annotation(tss_file="refseq.refflat.tss.bed")
        # Annotate peaks with genomic regions
        analysis.get_peak_genomic_location(genome=genome)
        # Annotate peaks with ChromHMM state
        if os.path.exists(os.path.join("data", "external", "HAP1_12_segments.annotated.bed")):
            analysis.get_peak_chromatin_state(chrom_state_file="data/external/HAP1_12_segments.annotated.bed")
        # Annotate peaks with closest gene, chromatin state,
        # genomic location, mean and variance measurements across samples
        analysis.annotate()
        attrs = (
            ['sample_name'] + sample_attributes
            if "sample_name" not in sample_attributes else sample_attributes)
        analysis.annotate_with_sample_metadata(attributes=attrs)
        analysis.to_pickle()

        # QC plots
        # plot general peak set features
        # analysis.plot_peak_characteristics()
        # plot coverage features across peaks/samples
        # analysis.plot_coverage()
        # analysis.plot_variance()

        quant_matrix = "accessibility"
        feature_name = "sites"

    if data_type == "RNA-seq":
        # Get gene expression
        analysis.get_gene_expression(
            samples=analysis.samples,
            sample_attributes=sample_attributes)

        quant_matrix = "expression_annotated"
        feature_name = "genes"

    # Unsupervised analysis
    unsupervised_analysis(
        analysis,
        quant_matrix=quant_matrix,
        attributes_to_plot=plotting_attributes,
        plot_max_attr=20,
        plot_max_pcs=6,
        plot_group_centroids=True,
        axis_ticklabels=False,
        axis_lines=True,
        always_legend=False,
        display_corr_values=False)

    # Supervised analysis
    comps = comparison_table[
        (comparison_table['data_type'] == data_type) &
        (comparison_table['comparison_type'] == 'differential')]
    if comps.empty:
        print("Comparison table has no comparisons with 'data_type'=='{}' and 'comparison_type'=='differential'.".format(data_type))
        print("Not performing differential analysis for this data type.")
        return 0

    if not (comps['sample_name'].value_counts() > 1).any():
        print("Detected simple design for differential analysis. Each sample belongs only to one group.")
        analysis.differential_results = differential_analysis(
            analysis,
            comps,
            data_type=data_type,
            samples=[s for s in analysis.samples if s.name in comps['sample_name'].tolist()],
            covariates=None,
            alpha=0.05,
            overwrite=True)
        analysis.differential_results = analysis.differential_results.set_index("index")
    else:
        print("Complex design for differential analysis. There are sample(s) belonging to more than one group.")
        print("Performing analysis independently for each comparison.")
        analysis.differential_results = pd.DataFrame()
        for comparison in comps['comparison_name'].unique():
            comp = comps[comps['comparison_name'] == comparison]
            res = differential_analysis(
                analysis,
                comp,
                data_type=data_type,
                samples=[s for s in analysis.samples if s.name in comp['sample_name'].tolist()],
                covariates=None, alpha=alpha, overwrite=True)
            analysis.differential_results = analysis.differential_results.append(res, ignore_index=True)
        analysis.differential_results = analysis.differential_results.set_index("index")
        analysis.differential_results.to_csv(
            os.path.join(analysis.results_dir, "differential_analysis_{}".format(data_type),
                         "differential_analysis.deseq_result.all_comparisons.csv"), index=True)
    analysis.to_pickle()

    diff = analysis.differential_results[
            (analysis.differential_results['padj'] < alpha) &
            (analysis.differential_results['log2FoldChange'].abs() > abs_fold_change)]
    if diff.empty:
        print("Differential analysis contains no significant {} at alpha {} and absolute fold change {}.".format(feature_name, alpha, abs_fold_change))
        return 0

    if diff.groupby('comparison_name').count().shape[0] > 1:
        differential_overlap(
            diff,
            getattr(analysis, quant_matrix).shape[0],
            data_type=data_type)

    plot_differential(
        analysis,
        analysis.differential_results,
        matrix=getattr(analysis, quant_matrix),
        comparison_table=comps,
        data_type=data_type,
        alpha=alpha,
        corrected_p_value=True,
        fold_change=abs_fold_change,
        rasterized=True,
        robust=True,
        group_wise_colours=True,
        group_variables=plotting_attributes)

    differential_enrichment(
        analysis,
        analysis.differential_results[
            (analysis.differential_results['padj'] < alpha) &
            (analysis.differential_results['log2FoldChange'].abs() > abs_fold_change)],
        data_type=data_type,
        genome=genome,
        directional=True,
        max_diff=1000,
        sort_var="pvalue",
        as_jobs=False)

    collect_differential_enrichment(
        analysis.differential_results[
            (analysis.differential_results['padj'] < alpha) &
            (analysis.differential_results['log2FoldChange'].abs() > abs_fold_change)],
        directional=True,
        data_type=data_type,
        permissive=False)

    if data_type == "RNA-seq":
        enrichment_table = pd.read_csv(
            os.path.join("{}/differential_analysis_{}".format(
                analysis.results_dir, data_type), "differential_analysis.enrichr.csv"))

        plot_differential_enrichment(
            enrichment_table,
            "enrichr",
            data_type=data_type,
            output_prefix="differential_analysis",
            direction_dependent=True,
            top_n=5)

    elif data_type in ["ATAC-seq", "ChIP-seq"]:
        for enrichment_name, enrichment_type in [
                ('motif', 'meme_ame'), ('lola', 'lola'), ('enrichr', 'enrichr')]:
            try:
                enrichment_table = pd.read_csv(
                    os.path.join("{}/differential_analysis_{}".format(
                        analysis.results_dir, data_type), "differential_analysis" + ".{}.csv".format(enrichment_type)))
            except pd.errors.EmptyDataError:
                print("Enrichment dataframe of {} is empty.".format(enrichment_type))
                continue

            plot_differential_enrichment(
                enrichment_table,
                enrichment_name,
                data_type=data_type,
                direction_dependent=True,
                top_n=5 if enrichment_name != "motif" else 300)

    return analysis


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
