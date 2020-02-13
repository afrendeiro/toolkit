#!/usr/bin/env python

"""
A helper script to run Location Overlap Analysis (LOLA)
of a single region set in various sets of region-based annotations.


Software requirements:

 * LOLA
"""

import os
import sys

from argparse import ArgumentParser
from ngs_toolkit.general import lola


def parse_arguments():
    """
    Argument Parsing.
    """
    parser = ArgumentParser(
        prog="python -m ngs_toolkit.recipes.lola", description=__doc__)
    parser.add_argument(dest="bed_file", help="BED file with query set regions.")
    parser.add_argument(
        dest="universe_file",
        help="BED file with universe where the query set came from.",
    )
    parser.add_argument(
        dest="output_folder", help="Output directory for produced files."
    )
    parser.add_argument(dest="genome", help="Genome assembly of the region set.")
    parser.add_argument(
        "--no-overwrite",
        action="store_false",
        help="Don't overwrite existing output files.",
    )
    parser.add_argument(
        "-c",
        "--cpus",
        dest="cpus",
        help="Number of CPUS/threads to use for analysis.",
        type=int,
    )
    return parser


def main(cli=None):
    print("LOLA analysis")
    args = parse_arguments().parse_args(cli)

    output_file = os.path.join(args.output_folder, "allEnrichments.tsv")
    if os.path.exists(output_file) and (not args.overwrite):
        print("Output exists and `overwrite` is False, so not doing anything.")
        return 0

    print("Starting LOLA analysis.")

    lola(
        args.bed_file,
        args.universe_file,
        args.output_folder,
        args.genome,
        cpus=args.cpus,
    )

    print("Done.")


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
