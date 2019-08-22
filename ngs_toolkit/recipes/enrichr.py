#!/usr/bin/env python

"""
A helper script to run enrichment analysis using
the Enrichr API on a gene set.
"""

import argparse
import os
import sys

import pandas as pd
from ngs_toolkit.general import enrichr


def parse_arguments():
    """
    Argument Parsing.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        dest="input_file",
        help="Input file with gene names.",
    )
    parser.add_argument(
        dest="output_file", help="Output CSV file with results."
    )
    parser.add_argument(
        "-a", "--max-attempts", type=int, default=5,
        dest="max_attempts", help="Maximum attempts to retry the API before giving up."
    )
    parser.add_argument(
        "--no-overwrite", action="store_false",
        dest="overwrite", help="Whether results should not be overwritten if existing."
    )
    args = parser.parse_args()

    return args


def main():
    print("Enrichr analysis")
    args = parse_arguments()
    if os.path.exists(args.output_file) and (not args.overwrite):
        print("Output exists and `overwrite` is False, so not doing anything.")
        return 0

    print("Reading input file.")

    df = pd.read_csv(args.input_file, header=None, names=["gene_name"])

    print("Starting Enrichr analysis.")
    res = enrichr(df, gene_set_libraries=None, kind="genes", max_attempts=args.max_attempts)

    print("Saving results.")
    res.to_csv(args.output_file, index=False)

    print("Done.")


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
