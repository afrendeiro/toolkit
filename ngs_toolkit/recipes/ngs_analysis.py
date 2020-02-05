#!/usr/bin/env python

"""
Perform full end-to-end analysis of ATAC-seq, ChIP-seq or RNA-seq data.

Produces quantification matrices, normalizes them,
performes unsupervised and supervised analysis as
well as enrichment analyisis of differential features,
all accompaigned with powerful visualizations.

Supervised analysis will only be performed if PEP configuration file contains a
`comparison table <https://ngs-toolkit.readthedocs.io/en/latest/comparison_table.html>`_ field.

In addition, this recipe uses variables provided in the project configuration
file ``project_name``, ``sample_attributes`` and ``group_attributes``.
"""


import os
import sys

from argparse import ArgumentParser

import matplotlib
import seaborn as sns

import peppy

from ngs_toolkit.atacseq import ATACSeqAnalysis
from ngs_toolkit.chipseq import ChIPSeqAnalysis
from ngs_toolkit.rnaseq import RNASeqAnalysis


# Set settings
sns.set(context="paper", style="ticks", palette="colorblind", color_codes=True)
matplotlib.rc("text", usetex=False)


def parse_arguments():
    """
    Global options for analysis.
    """
    parser = ArgumentParser(
        prog="python -m ngs_toolkit.recipes.ngs_analysis", description=__doc__)
    parser.add_argument(
        dest="config_file", help="YAML project configuration file.", type=str)
    parser.add_argument(
        "-n",
        "--analysis-name",
        dest="name",
        default=None,
        help="Name of analysis. Will be the prefix of output_files. "
        "By default it will be the name of the Project given in the YAML configuration.",
        type=str,
    )
    parser.add_argument(
        "-o",
        "--results-output",
        default="results",
        dest="results_dir",
        help="Directory for analysis output files. "
             "Default is 'results' under the project roort directory.",
        type=str,
    )
    parser.add_argument(
        "-t",
        "--data-type",
        default=None,
        choices=["ATAC-seq", "RNA-seq", "ChIP-seq"],
        dest="data_type",
        help="Data type to restrict analysis to. "
             "Default is to run separate analysis for each data type.",
        type=str,
    )
    parser.add_argument(
        "-q",
        "--pass-qc",
        action="store_true",
        dest="pass_qc",
        help="Whether only samples with a 'pass_qc' value of '1' "
             "in the annotation sheet should be used.",
    )
    parser.add_argument(
        "-a", "--alpha", default=0.05, dest="alpha",
        help="Alpha value of confidence for supervised analysis.", type=str
    )
    parser.add_argument(
        "-f",
        "--fold-change",
        default=0,
        dest="abs_fold_change",
        help="Absolute log2 fold change value for supervised analysis.",
        type=str,
    )
    return parser


def main(cli=None):
    args = parse_arguments().parse_args(cli)

    # Start project
    print("Starting peppy project with project"
          "configuration file: '{}'".format(args.config_file))
    prj = peppy.Project(args.config_file)
    print("Changing directory to project root"
          "directory: '{}'.".format(prj.metadata.output_dir))
    os.chdir(prj.metadata.output_dir)
    if args.pass_qc:
        print("Filtering samples out which didn't pass QC"
              "as specified in sample annotation in column 'pass_qc'")
        prj._samples = [
            s for s in prj._samples
            if s.pass_qc not in ["0", 0, "False", False]]

    # ANALYSIS
    if args.data_type is None:
        print(
            "Type of analysis not specified. Will run independent analysis"
            "for all types of data in the sample annotation sheet."
        )
        data_types = sorted(list(set([s.protocol for s in prj._samples])))
        print("Sample data types: '{}'.".format(",".join(data_types)))
    else:
        print("Type of analysis specified. Will run only"
              "analysis for samples of type '{}'.".format(args.data_type))
        data_types = [args.data_type]
        print("Sample data types: '{}'.".format(",".join(data_types)))
    if args.name is None:
        print(
            "Analysis name not specified, will use name in"
            "project configuration file: '{}'.".format(prj.project_name)
        )
        args.name = prj.project_name

    for data_type in data_types:
        print("Starting analysis for samples of type: '{}'.".format(data_type))
        samples = [s for s in prj._samples if (s.protocol == data_type)]
        if len(samples) > 0:
            print(
                "Samples under consideration: '{}'. ".format(",".join([s.name for s in samples]))
                + "Total of {} samples.".format(len([s.name for s in samples]))
            )
        else:
            raise ValueError("There were no valid samples for this analysis type!")

        kwargs = {"prj": prj, "samples": samples, "results_dir": args.results_dir}
        if data_type in ["ATAC-seq"]:
            print("Initializing ATAC-seq analysis")
            analysis = ATACSeqAnalysis(
                name=args.name + "_atacseq", **kwargs
            )
        elif data_type in ["ChIP-seq"]:
            print("Initializing ChIP-seq analysis")
            analysis = ChIPSeqAnalysis(
                name=args.name + "_chipseq", **kwargs
            )
        elif data_type in ["RNA-seq"]:
            print("Initializing RNA-seq analysis")
            analysis = RNASeqAnalysis(
                name=args.name + "_rnaseq", **kwargs
            )

        print("Running main analysis.")
        main_analysis_pipeline(
            analysis, alpha=args.alpha, abs_fold_change=args.abs_fold_change)
        print("`ngs_analysis` recipe completed successfully!")


def main_analysis_pipeline(a, alpha=0.05, abs_fold_change=0):
    # TODO: annotate with chromatin state
    # TODO: handle the genome vs transcriptome ambiguity

    genomes = list(set(s.genome for s in a.samples))

    if len(genomes) != 1:
        raise ValueError(
            "Samples under analysis have more than"
            "one genome assembly: '{}'.".format("', '".join(genomes))
        )

    if isinstance(a, ATACSeqAnalysis):

        # GET CONSENSUS PEAK SET, ANNOTATE IT, PLOT
        # Get consensus peak set from all samples
        a.get_consensus_sites()
        a.calculate_peak_support()

        # GET CHROMATIN OPENNESS MEASUREMENTS, PLOT
        # Get coverage values for each peak in each sample
        a.measure_coverage()
        # normalize coverage values
        a.normalize(method="vst")

        # Annotate peaks with closest gene
        a.get_peak_gene_annotation()
        # Annotate peaks with genomic regions
        a.get_peak_genomic_location()
        # Annotate peaks with chromatin state

    if isinstance(a, RNASeqAnalysis):
        # Get gene expression
        a.get_gene_expression()

    # Annotate peaks with closest gene, chromatin state,
    # genomic location, mean and variance measurements across samples
    a.annotate_features()
    a.to_pickle()

    # Unsupervised analysis
    a.unsupervised_analysis(
        plot_max_attr=20,
        plot_max_pcs=6,
        plot_group_centroids=True,
        axis_ticklabels=False,
        axis_lines=True,
        always_legend=False,
        display_corr_values=False,
    )

    # Supervised analysis
    if a.comparison_table.empty:
        print(
            "Comparison table has no comparisons with 'data_type'=='{}'"
            "and 'comparison_type'=='differential'.".format(
                a.data_type
            )
        )
        print("Not performing differential analysis for this data type.")
        a.generate_report(pip_versions=True)
        a.to_pickle()
        return

    a.differential_analysis()
    a.to_pickle()

    diff = a.differential_results[
        (a.differential_results["padj"] < alpha)
        & (a.differential_results["log2FoldChange"].abs() > abs_fold_change)
    ]
    if diff.empty:
        print(
            "Differential analysis contains no significant {}s"
            "at alpha {} and absolute fold change {}.".format(
                a.var_unit_name, alpha, abs_fold_change
            )
        )
        a.generate_report(pip_versions=True)
        a.to_pickle()
        return

    if diff.groupby("comparison_name").count().shape[0] > 1:
        a.differential_overlap(diff)

    a.plot_differential(
        alpha=alpha,
        corrected_p_value=True,
        fold_change=abs_fold_change,
        rasterized=True,
        robust=True,
        group_wise_colours=True,
    )

    a.differential_enrichment(
        # TODO: have a way to automatically check what is callable
        steps=['enrichr'],
        directional=True,
        max_diff=1000,
        sort_var="pvalue",
        distributed=False)

    # TODO: is this actually needed? vvv
    # a.collect_differential_enrichment(directional=True, permissive=False)

    a.plot_differential_enrichment(direction_dependent=True, top_n=5)

    a.generate_report(pip_versions=True)
    a.to_pickle()


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
