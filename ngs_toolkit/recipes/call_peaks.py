#!/usr/bin/env python

"""
This is the "call_peaks" recipe for ngs_toolkit.

This will call peaks for ChIP-seq samples given a comparison table mapping relationships between samples.
"""


from argparse import ArgumentParser
import os
import sys

import pandas as pd

from peppy import Project
from ngs_toolkit.chipseq import ChIPSeqAnalysis


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
        "-c", "--comparison-table",
        dest="comparison_table",
        default=None,
        help="Comparison table to use for peak calling. If not provided will use a file"
             "named `comparison_table.csv` in the same directory of the given YAML Project configuration file.",
        type=str)
    parser.add_argument(
        "-t", "--only-toggle",
        action="store_true",
        dest="only_toggle",
        help="Whether only comparisons with 'toggle' value of '1' "
        "in the should be performed.")
    parser.add_argument(
        "-qc", "--pass-qc",
        action="store_true",
        dest="pass_qc",
        help="Whether only samples with a 'pass_qc' attribute should be included."
        " Default is False.")
    parser.add_argument(
        "-j", "--as-jobs",
        action="store_true",
        dest="as_job",
        help="Whether jobs should be created for each sample, or "
        "it should run in serial mode.")
    parser.add_argument(
        "-o", "--results-output",
        default="results",
        dest="results_dir",
        help="Directory for analysis output files. "
        "Default is 'results' under the project root directory.",
        type=str)
    return parser


def main():
    parser = ArgumentParser(
        prog="call_peaks_recipe",
        description="Call peaks recipe."
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

    # ANALYSIS
    data_types = sorted(list(set([s.library for s in prj.samples])))
    print("Sample data types: '{}'.".format(",".join(data_types)))

    if args.comparison_table is None:
        print("Comparison table not specified, will use name in project configuration file: '{}'.".format(prj.project_name))
        args.comparison_table = os.path.join(os.path.dirname(args.config_file), "comparison_table.csv")

    # Read comparison table
    try:
        comparison_table = pd.read_csv(args.comparison_table)
    except IOError as e:
        print("Comparison table could not be opened: {}".format(args.comparison_table))
        raise e

    comparison_table = comparison_table[comparison_table['comparison_type'] == 'peaks']

    if args.only_toggle:
        print("Filtering out comparisons marked with toggle != 1")
        comparison_table = comparison_table[comparison_table['toggle'] == 1]

    comps = comparison_table["comparison_name"].unique()
    if len(comps) > 0:
        print(
            "comparisons under consideration: '{}'. ".format(",".join(comps)) +
            "Total of {} comparisons.".format(len(comps)))
    else:
        raise ValueError("There were no valid comparisons in the comparison table!")

    for data_type in [dt for dt in data_types if dt in ["ChIP-seq", "ChIPmentation"]]:
        print("Starting analysis for samples of type: '{}'.".format(data_type))
        samples = [s for s in prj.samples if (s.library == data_type)]
        if len(samples) > 0:
            print(
                "Samples under consideration: '{}'. ".format(",".join([s.name for s in samples])) +
                "Total of {} samples.".format(len([s.name for s in samples])))
        else:
            raise ValueError("There were no valid samples for this analysis type!")

        print("Initializing ChIP-seq analysis")
        if hasattr(prj, "project_name"):
            name = prj.project_name
        else:
            name = os.path.basename(os.path.abspath(os.curdir))
        analysis = ChIPSeqAnalysis(
            name=name + "_chipseq", prj=prj,
            samples=samples, results_dir=args.results_dir)

        # Call peaks
        analysis.call_peaks_from_comparisons(comparison_table, as_jobs=args.as_jobs)

        # # Get summary of peak calls
        # peak_counts = analysis.summarize_peaks_from_comparisons(comparison_table)
        # peak_counts.to_csv(os.path.join("results_pipeline", "chipseq_peaks", "peak_count_summary.csv"), index=False)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
