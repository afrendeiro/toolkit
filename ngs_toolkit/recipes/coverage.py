#!/usr/bin/env python

"""
A helper script to calculate the read coverage of a BAM file
in regions from a BED file.
Ensures the same order and number of lines as input BED file.

Software requirements:

 * None
"""

import os
import sys

from argparse import ArgumentParser

import pandas as pd

from ngs_toolkit.utils import count_reads_in_intervals
from ngs_toolkit.utils import read_bed_file_three_columns
from ngs_toolkit.utils import to_bed_index


def parse_arguments():
    """
    Argument Parsing.
    """
    parser = ArgumentParser(
        prog="python -m ngs_toolkit.recipes.coverage", description=__doc__)
    parser.add_argument(
        dest="bed_file",
        help="Input BED file with regions to quantify.",
    )
    parser.add_argument(
        dest="bam_file",
        help="Input BAM file with reads.",
    )
    parser.add_argument(
        dest="output_bed", help="Output BED file with counts for each region."
    )
    parser.add_argument(
        "--no-overwrite", action="store_false",
        dest="overwrite",
        help="Whether results should not be overwritten if existing."
    )
    return parser


def main(cli=None):
    """Measure coverage of BAM file in BED file regions."""
    print("Parsing CLI.")
    args = parse_arguments().parse_args(cli)

    if os.path.exists(args.output_bed) and (not args.overwrite):
        print("Output exists and `overwrite` is False, so not doing anything.")
        return 0

    print("Getting regions.")
    sites_str = to_bed_index(args.bed_file)
    print("Quantifying.")
    res = count_reads_in_intervals(args.bam_file, sites_str)

    print("Merging with input set.")
    # make sure there is an entry for each region in input file
    input_bed = read_bed_file_three_columns(args.bed_file).set_index("name")
    res = input_bed.join(pd.Series(res, name="sample")).fillna(0)
    res.loc[:, "sample"] = res.loc[:, "sample"].astype(int)

    print("Saving results.")
    res.to_csv(args.output_bed, index=False, header=False, sep="\t")

    print("Done.")


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
