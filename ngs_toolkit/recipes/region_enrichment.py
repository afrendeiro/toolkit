#!/usr/bin/env python

"""
A helper script to run enrichment analysis
of a single region set in region-based set of annotations.
"""

import os
import pandas as pd
import sys

from argparse import ArgumentParser
from ngs_toolkit.atacseq import ATACSeqAnalysis
from ngs_toolkit.utils import bed_to_index


def parse_arguments():
    """
    Argument Parsing.
    """
    parser = ArgumentParser(
        prog="python -m ngs_toolkit.recipes.region_enrichment", description=__doc__)
    parser.add_argument(
        dest="bed_file",
        help="BED file with regions.")
    parser.add_argument(
        dest="pep",
        help="The analysis' PEP config file.")
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
    return parser


def main(cli=None):
    print("Region type analysis")
    # Parse command-line arguments.
    args = parse_arguments().parse_args(cli)
    if os.path.exists(args.output_file) and (not args.overwrite):
        print("Output exists and `overwrite` is False, so not doing anything.")
        return 0

    print("Reading up the analysis object.")
    a = ATACSeqAnalysis(from_pep=args.pep)
    a.load_data()
    # (
    #     "genomic_region",
    #     "region_annotation_mapping",
    #     "region_annotation_b_mapping",
    # ),
    # (
    #     "chromatin_state",
    #     "chrom_state_annotation_mapping",
    #     "chrom_state_annotation_b_mapping",
    # ),
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
