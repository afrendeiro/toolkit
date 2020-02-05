#!/usr/bin/env python

"""
Call peaks for ChIP-seq samples given a comparison table
mapping foreground-background relationships between samples.
"""


import sys

from argparse import ArgumentParser

import pandas as pd

from ngs_toolkit.chipseq import ChIPSeqAnalysis


def parse_arguments():
    """
    Global options for analysis.
    """
    parser = ArgumentParser(
        prog="python -m ngs_toolkit.recipes.call_peaks", description=__doc__)
    parser.add_argument(
        dest="config_file", help="YAML project configuration file.", type=str
    )
    parser.add_argument(
        "-c",
        "--comparison-table",
        dest="comparison_table",
        default=None,
        help="Comparison table to use for peak calling. If not provided will use a file"
        "named `comparison_table.csv` in the same directory of the given YAML Project configuration file.",
        type=str,
    )
    parser.add_argument(
        "-t",
        "--only-toggle",
        action="store_true",
        dest="only_toggle",
        help="Whether only comparisons with 'toggle' value of '1' or 'True' should be performed.",
    )
    parser.add_argument(
        "-qc",
        "--pass-qc",
        action="store_true",
        dest="pass_qc",
        help="Whether only samples with a 'pass_qc' attribute should be included."
        " Default is :obj:`False`.",
    )
    parser.add_argument(
        "-j",
        "--as-jobs",
        action="store_true",
        dest="as_job",
        help="Whether jobs should be created for each sample, or "
        "it should run in serial mode.",
    )
    parser.add_argument(
        "-o",
        "--results-output",
        default="results",
        dest="results_dir",
        help="Directory for analysis output files. "
        "Default is 'results' under the project root directory.",
        type=str,
    )
    return parser


def main(cli=None):
    args = parse_arguments().parse_args(cli)

    # Analysis
    print(
        "Starting Analysis from PEP configuration file: '{}'".format(args.config_file)
    )
    analysis = ChIPSeqAnalysis(
        from_pep=args.config_file, results_dir=args.results_dir
    )
    chip_data_types = ["ChIP-seq", "ChIPmentation"]
    analysis.samples = [s for s in analysis.samples if s.protocol == chip_data_types]

    # Samples
    # # filter QC if needed
    if args.pass_qc:
        analysis.samples = [
            s for s in analysis.samples if s.pass_qc not in ["0", 0, "False", False]
        ]
    if analysis.samples:
        print(
            "Samples under consideration: '{}'. ".format(
                ",".join([s.name for s in analysis.samples])
            )
            + "Total of {} samples.".format(len([s.name for s in analysis.samples]))
        )
    else:
        raise ValueError("There were no valid samples for this analysis type!")

    # Comparison table
    # # add provided
    if args.comparison_table is not None:
        analysis.comparison_table = pd.read_csv(args.comparison_table)
    # # or make sure analysis has one
    else:
        if not hasattr(analysis, "comparison_table"):
            raise ValueError(
                "Analysis doesn't have a 'comparison_table' and this was not provided."
            )

    # # filter comparisons if needed
    if args.only_toggle:
        print("Filtering out comparisons marked with toggle != 1")
        analysis.comparison_table = analysis.comparison_table[
            analysis.comparison_table["toggle"] == 1
        ]

    comps = analysis.comparison_table["comparison_name"].unique()
    if comps:
        print(
            "comparisons under consideration: '{}'. ".format(",".join(comps))
            + "Total of {} comparisons.".format(len(comps))
        )
    else:
        raise ValueError("There were no valid comparisons in the comparison table!")

    # Call peaks
    analysis.call_peaks_from_comparisons(distributed=args.as_jobs)

    # # Get summary of peak calls
    # peak_counts = analysis.summarize_peaks_from_comparisons(comparison_table)
    # peak_counts.to_csv(os.path.join("results_pipeline", "chipseq_peaks", "peak_count_summary.csv"), index=False)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
