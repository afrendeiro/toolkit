#!/usr/bin/env python

"""
A helper script to calculate the read coverage of a BAM file
in regions from a BED file.
Ensures the same order and number of lines as input BED file.
"""

import argparse
import os
import sys

import pandas as pd
from ngs_toolkit.utils import to_bed_index, count_reads_in_intervals, \
    read_bed_file_three_columns


def parse_arguments():
    """
    Argument Parsing.
    """
    parser = argparse.ArgumentParser()
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
    args = parser.parse_args()

    return args


def main() -> int:
    """Measure coverage of BAM file in BED file regions."""
    print("Enrichr analysis")
    args = parse_arguments()
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

    print("Saving results.")
    res.to_csv(args.output_bed, index=False, header=False, sep="\t")

    print("Done.")


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
