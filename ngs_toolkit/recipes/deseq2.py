#!/usr/bin/env python

"""
Perform differential expression using DESeq2
by comparing sample groups using a formula.

Software requirements:

 * DESeq2
"""

import os
import pandas as pd
import sys

from argparse import ArgumentParser
from ngs_toolkit.general import deseq_analysis


def parse_arguments():
    parser = ArgumentParser(
        prog="python -m ngs_toolkit.recipes.deseq2", description=__doc__)
    parser.add_argument(
        dest="work_dir",
        help="Working directory. Should contain required files for DESeq2.")
    parser.add_argument(
        "--output-prefix",
        dest="output_prefix",
        default="differential_analysis",
        type=str,
        help="Prefix for output files.")
    parser.add_argument(
        "--formula",
        default="~ sample_group",
        type=str,
        help="R-style formula for differential expression. Defaults to '~ sample_group'.")
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
    parser.add_argument(
        "--no-save-inputs",
        action="store_false",
        default=True,
        help="Don't write inputs to disk.")

    # To enable the loop to pass args directly on to the pipelines...

    # args = parser.parse_args("--output_prefix differential_analysis --formula '~sample_group' --overwrite /scratch/lab_bock/shared/projects/baf-time_course/results/differential_analysis_ATAC-seq/ARID2_KO".split(" "))
    return parser


def main(cli=None):
    args = parse_arguments().parse_args(cli)

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
        create_subdirectories=False,
        save_inputs=not args.no_save_inputs)

    print("Found {} differentially expressed genes with p < {}.".format(
        res[res['pvalue'] < args.alpha].shape[0], args.alpha))
    print("Found {} differentially expressed genes with FDR < {}.".format(
        res[res['padj'] < args.alpha].shape[0], args.alpha))


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
