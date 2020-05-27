#!/usr/bin/env python

"""
Compute fraction of reads in peaks (FRiP) based on a consensus set of regions
derived from several samples.

A consensus region set can be passed, otherwise it will either try to use an
existing one for that analysis or produce one on the fly.


Software requirements:

 * awk
 * samtools
"""


import os
import sys

from argparse import ArgumentParser

import pybedtools

from ngs_toolkit.atacseq import ATACSeqAnalysis
from ngs_toolkit.chipseq import ChIPSeqAnalysis


def parse_arguments():
    """
    Global options for analysis.
    """
    parser = ArgumentParser(
        prog="python -m ngs_toolkit.recipes.region_set_frip", description=__doc__
    )
    parser.add_argument(dest="config_file", help="YAML project configuration file.", type=str)
    parser.add_argument(
        "-r",
        "--region-set",
        dest="region_set",
        default=None,
        help="BED file with region set derived from several samples or Oracle region set. "
        "If unset, will try to get the `sites` attribute of an existing analysis object "
        "if existing, otherwise will create a region set from the peaks of all samples.",
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
        "--computing-configuration",
        dest="computing_configuration",
        help="Which `divvy` computing configuration to use for distributed jobs."
        " Type divvy list to see all options. Defaults to the value in the "
        "ngs_toolkit configuration.",
    )
    parser.add_argument(
        "--permissive",
        action="store_true",
        dest="permissive",
        help="If creating regions set, allow sample files to be missing and use what is present.",
    )
    return parser


def main(cli=None):
    args = parse_arguments().parse_args(cli)

    for data_type, clax in [
        ("ATAC-seq", ATACSeqAnalysis),
        ("ChIP-seq", ChIPSeqAnalysis),
    ]:
        an = clax(from_pep=args.config_file)

        if not an.samples:
            continue

        if args.pass_qc:
            an.samples = [s for s in an.samples if getattr(s, "pass_qc", None) in ["1", "1.0", 1]]

        if data_type == "ChIP-seq" and not hasattr(an, "comparison_table"):
            msg = (
                "ChIP-seq analysis must have comparison_table specified in "
                "the project config in order to relate"
                " foreground and backgound sample groups."
            )
            print(msg)
            raise ValueError(msg)

        if args.region_set is not None:
            print("Loading given region set: '{}'".format(args.region_set))
            an.sites = pybedtools.BedTool(args.region_set)
        else:
            print("Trying to load existing consensus region set.")
            an.load_data(only_these_keys=["sites"])

        if not hasattr(an, "sites"):
            print("Not found. Producing a new consensus region set.")
            an.get_consensus_sites(permissive=args.permissive)
        else:
            print("Using region set in BED format: '{}'".format(an.sites.fn))

        calculate_region_set_frip(
            region_set=an.sites.fn,
            samples=an.samples,
            computing_configuration=args.computing_configuration,
        )


def calculate_region_set_frip(region_set, samples, computing_configuration=None):
    """
    """
    from ngs_toolkit.utils import submit_job

    for sample in samples:
        sample.sample_root = os.path.join(
            sample.project.root_dir, sample.project._config.results_subdir, sample.name
        )
        inside_reads = os.path.join(sample.sample_root, "region_set_frip.inside_reads.txt")
        all_reads = os.path.join(sample.sample_root, "region_set_frip.all_reads.txt")

        job_name = sample.name + ".region_set_frip"
        log_file = os.path.join(sample.sample_root, job_name + ".log")
        job_file = os.path.join(sample.sample_root, job_name + ".sh")
        sample_stats = os.path.join(sample.sample_root, "stats.tsv")

        cmd = "\n".join(
            [
                """samtools view -c -L {} {} > {}""".format(
                    region_set, sample.aligned_filtered_bam, inside_reads
                ),
                """samtools view -c {} > {}""".format(sample.aligned_filtered_bam, all_reads),
                'calc(){ awk "BEGIN { print "$*" }"; }',
                "IN=`cat {}`".format(inside_reads),
                "ALL=`cat {}`".format(all_reads),
                "FRIP=`calc $IN/$ALL`",
                'echo "region_set_frip\\t$FRIP\\t." >> {}'.format(sample_stats),
                "date",
            ]
        )
        submit_job(
            cmd,
            job_file,
            log_file,
            jobname=job_name,
            computing_configuration=computing_configuration,
        )


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
