#!/usr/bin/env python

"""
A helper script to run enrichment analysis
of a single region set in region-based set of annotations.
"""

import argparse
import os
import sys
from ngs_toolkit.atacseq import ATACSeqAnalysis
from ngs_toolkit.utils import bed_to_index
import pandas as pd


def parse_arguments():
    """
    Argument Parsing.
    """
    parser = argparse.ArgumentParser(
        description=__doc__)
    parser.add_argument(
        dest="bed_file",
        help="BED file with regions.")
    parser.add_argument(
        dest="pickle",
        help="The analysis' pickle object.")
    parser.add_argument(
        "--output-file",
        dest="output_file",
        default="region_type_enrichment.csv",
        type=str,
        help="Output file.")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help="Don't overwrite any existing directory or file.")
    args = parser.parse_args()

    return args


def main():
    print("Region type analysis")
    # Parse command-line arguments.
    args = parse_arguments()
    if os.path.exists(args.output_file) and (not args.overwrite):
        print("Output exists and `overwrite` is False, so not doing anything.")
        return 0

    print("Reading up the analysis object.")
    a = ATACSeqAnalysis(pickle_file=args.pickle, from_pickle=True)
    print("Reading up the BED file.")
    df = pd.read_csv(args.bed_file, sep="\t", header=None)
    df.columns = ['chrom', 'start', 'end']
    print("Getting the index.")
    index = bed_to_index(df)
    print("Doing enrichment.")
    enr = a.region_context_enrichment(index)
    print("Saving.")
    enr.to_csv(args.output_file)
    print("Done.")


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
