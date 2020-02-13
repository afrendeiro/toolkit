#!/usr/bin/env python

"""
A helper script to run enrichment analysis using the Enrichr API on a gene set.

Software requirements: None
"""

from argparse import ArgumentParser
import os
import sys

from ngs_toolkit.general import enrichr


def parse_arguments():
    """
    Argument Parsing.
    """
    parser = ArgumentParser(
        prog="python -m ngs_toolkit.recipes.enrichr", description=__doc__
    )
    parser.add_argument(
        dest="input_file",
        help="Input file with a gene name per row and no header.",
    )
    parser.add_argument(
        dest="output_file", help="Output CSV file with results."
    )
    parser.add_argument(
        "-a",
        "--max-attempts",
        type=int,
        default=5,
        dest="max_attempts",
        help="Maximum attempts to retry the API before giving up.",
    )
    parser.add_argument(
        "--no-overwrite",
        action="store_false",
        dest="overwrite",
        help="Whether results should not be overwritten if existing.",
    )
    return parser


def main(cli=None):
    print("Enrichr analysis")
    args = parse_arguments().parse_args(cli)

    if os.path.exists(args.output_file) and (not args.overwrite):
        print("Output exists and `overwrite` is False, so not doing anything.")
        return 0

    print("Reading input file.")

    with open(args.input_file, "r") as handle:
        genes = handle.readlines()

    print("Found {} genes in input.".format(len(genes)))

    print("Starting Enrichr analysis.")
    res = enrichr(
        genes,
        gene_set_libraries=None,
        kind="genes",
        max_attempts=args.max_attempts,
    )

    print("Saving results.")
    res.to_csv(args.output_file, index=False)

    print("Done.")


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
