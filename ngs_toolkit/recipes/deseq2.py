#!/usr/bin/env python

"""
Perform differential expression using DESeq2
by comparing sample groups using a formula.
"""

import argparse
import os
import sys
from ngs_toolkit.general import deseq_analysis
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=__doc__)
    parser.add_argument(
        dest="work_dir",
        help="Working directory. Should contain required files for DESeq2.")
    parser.add_argument(
        "--output_prefix",
        dest="output_prefix",
        default="differential_analysis",
        type=str,
        help="Prefix for output files.")
    parser.add_argument(
        "--formula",
        default="~ sample_group",
        type=str,
        help="R-style formula for differential expression. Default = '~ sample_group'.")
    parser.add_argument(
        "--alpha",
        default=0.05,
        type=float,
        help="Significance level to call differential expression. All results will be output anyway.")
    parser.add_argument(
        "-d",
        "--dry-run",
        action="store_true",
        help="Don't actually do anything.")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help="Don't overwrite any existing directory or file.")

    # To enable the loop to pass args directly on to the pipelines...
    args = parser.parse_args()

    return args


def main():
    args = parse_arguments()

    # sample annotation
    print("Reading experiment_matrix")
    experiment_matrix = pd.read_csv(
        os.path.join(args.work_dir, "experiment_matrix.csv"))
    # comparison table
    print("Reading comparison_matrix")
    comparison_table = pd.read_csv(
        os.path.join(args.work_dir, "comparison_table.csv"))
    # count matrix
    print("Reading count_matrix")
    count_matrix = pd.read_csv(
        os.path.join(args.work_dir, "count_matrix.csv"), index_col=0)

    print("Differential expression with DESeq2")
    res = deseq_analysis(
        count_matrix,
        experiment_matrix,
        comparison_table,
        formula=args.formula,
        output_dir=args.work_dir,
        output_prefix=args.output_prefix,
        overwrite=args.overwrite, alpha=args.alpha,
        create_subdirectories=False)

    print("Found {} differentially expressed genes with p < {}.\n".format(
        res[res['pvalue'] < args.alpha].shape[0], args.alpha))
    print("Found {} differentially expressed genes with FDR < {}.".format(
        res[res['padj'] < args.alpha].shape[0], args.alpha))


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
