#!/usr/bin/env python

"""
Merge signal from various ATAC-seq or ChIP-seq samples
given a set of attributes to group samples by.

It produces merged BAM and bigWig files for all signal in the samples but
is also capable of producing this for nucleosomal/nucleosomal free signal
based on fragment length distribution if data is paired-end sequenced.
This signal may optionally be normalized for each group.
It is also capable of parallelizing work in jobs.

Software requirements:

 * samtools
 * sambamba
 * deeptools
"""


import os
import sys

from argparse import ArgumentParser

import pandas as pd

from ngs_toolkit import _LOGGER
from ngs_toolkit import __version__
from ngs_toolkit import Analysis


def parse_arguments():
    """
    Global options for recipe.
    """
    parser = ArgumentParser(
        prog="python -m ngs_toolkit.recipes.merge_signal", description=__doc__)
    parser.add_argument(
        dest="config_file", help="YAML project configuration file.", type=str
    )
    parser.add_argument(
        "-a",
        "--attributes",
        dest="attributes",
        default=None,
        help="Attributes to merge samples by. A comma-delimited string with no spaces. "
        "By default will use values in the project config `group_attributes`.",
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
        "-j",
        "--as-jobs",
        action="store_true",
        dest="as_job",
        help="Whether jobs should be created for each sample, or "
        "it should run in serial mode.",
    )
    parser.add_argument(
        "--cpus",
        dest="cpus",
        default=8,
        help="CPUs/Threads to use per job if `--as-jobs` is on.",
    )
    parser.add_argument(
        "--normalization-method",
        dest="normalization_method",
        default="RPGC",
        help="Method to normalize tracks regarding sequenced depth. "
        "One of the methods in https://deeptools.readthedocs.io/en/develop/"
        "content/tools/bamCoverage.html#"
        "Read%20coverage%20normalization%20options",
    )
    parser.add_argument(
        "--nucleosome",
        action="store_true",
        dest="nucleosome",
        help="Whether to produce nucleosome/nucleosome-free signal files.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        dest="overwrite",
        help="Whether to overwrite existing files.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default="merged",
        dest="output_dir",
        help="Directory for output files. "
        "Default is 'merged' under the project root directory.",
        type=str,
    )
    parser.add_argument(
        "-d",
        "--dry-run",
        action="store_true",
        dest="dry_run",
        help="Whether to do everything except running commands.",
    )
    return parser


def main(cli=None):
    args = parse_arguments().parse_args(cli)
    _LOGGER.info(
        "This is the 'merge_signal' recipe from ngs_toolkit, "
        "version: %s", __version__
    )
    # Start project
    _LOGGER.debug(
        "Starting Analysis with PEP project configuration file: "
        "'%s'", args.config_file
    )
    an = Analysis(from_pep=args.config_file)
    if args.pass_qc:
        _LOGGER.info(
            "Filtering samples out which didn't pass QC as specified in sample "
            "annotation in column 'pass_qc'"
        )
        an.samples = [
            s for s in an.samples if getattr(s, "pass_qc") not in ["0", 0, "False", False]
        ]

    if an.samples:
        print(
            "Samples under consideration: '{}'. ".format(
                ",".join([s.name for s in an.samples])
            )
            + "\nTotal of {} samples.".format(len([s.name for s in an.samples]))
        )
    else:
        raise ValueError("There were no valid samples after filtering for quality!")

    # Get only samples with signal
    an.samples = [
        s for s in an.samples
        if getattr(s, "protocol", None)
        in ["ATAC-seq", "ChIP-seq", "ChIPmentation"]]
    an.set_organism_genome()
    sheet = an.prj.sheet.reindex([s.name for s in an.samples])

    _LOGGER.info(
        "Selecting samples with appropriate data type."
        "\nSamples under consideration: '%s'. "
        "\nTotal of %i samples.",
        ",".join(sheet["sample_name"].tolist()),
        sheet.shape[0]
    )

    # Get default attributes if not set
    if args.attributes is None:
        if an.group_attributes:
            args.attributes = an.group_attributes
        else:
            _LOGGER.error(
                "Sample attributes to group by were not set and none could be"
                " found in project configuration file!"
                " Aborting!"
            )
            return 1
    else:
        if "," in args.attributes:
            args.attributes = args.attributes.split(",")
        else:
            args.attributes = [args.attributes]

    _LOGGER.info(
        "Using the following attributes to merge samples: '%s', "
        "resulting in a total of %i groups.",
        "', '".join(args.attributes),
        len(sheet.groupby(args.attributes).groups.items())
    )

    merge_signal(
        sheet,
        an.samples,
        args.attributes,
        output_dir=args.output_dir,
        normalization_method=args.normalization_method,
        nucleosome=args.nucleosome,
        overwrite=args.overwrite,
        cpus=args.cpus,
        as_job=args.as_job,
        dry_run=args.dry_run,
    )


def merge_signal(
        sheet,
        samples,
        attributes,
        output_dir="merged",
        normalization_method="RPGC",
        nucleosome=False,
        overwrite=False,
        cpus=8,
        as_job=False,
        dry_run=False,
):
    """
    Merge signal for ``samples`` aggregated by the attributes in ``attributes``.
    """
    import re
    import subprocess

    def add_cmd(cmd, target, overwrite):
        if not overwrite:
            if not os.path.exists(target):
                return cmd
            else:
                return ""
        else:
            return cmd

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for attrs, index in sheet.groupby(attributes).groups.items():
        name = "_".join([a for a in attrs if not pd.isnull(a)])
        print(name)
        name = re.sub(r"_{2}", "_", name)
        name = re.sub(r"_$", "", name)
        sample_subset = [
            s for s in samples if s.name in sheet.loc[index, "sample_name"].tolist()
        ]
        bams = [s.aligned_filtered_bam for s in sample_subset]

        genomes = list(set(s.genome for s in sample_subset))
        if len(genomes) != 1:
            raise ValueError(
                "Samples under same attribute group have more than one genome assembly: '{}'.".format(
                    "', '".join(genomes)
                )
            )
        else:
            genome = genomes[0]

        job_file = os.path.join(output_dir, name + ".merge_signal.sh")
        log_file = os.path.join(output_dir, name + ".merge_signal.log")
        output_bam = os.path.join(output_dir, name + ".merged.bam")
        output_sorted_bam = os.path.join(output_dir, name + ".merged.sorted.bam")
        output_bigwig = os.path.join(output_dir, name + ".bigWig")
        output_nucleosome_free_reads = os.path.join(
            output_dir, name + ".nucleosome_free_reads.bam"
        )
        output_nucleosome_free_reads_bigwig = os.path.join(
            output_dir, name + ".nucleosome_free_reads.bigWig"
        )
        output_nucleosome_reads = os.path.join(
            output_dir, name + ".nucleosome_reads.bam"
        )
        output_nucleosome_reads_bigwig = os.path.join(
            output_dir, name + ".nucleosome_reads.bigWig"
        )

        # Get region databases from config
        _LOGGER.debug(
            "Getting genome assembly chromosome sizes for genome "
            "'%s' from configuration.",
            genome
        )

        # Prepare commands
        cmds = ["#!/bin/sh", "date"]

        # # merge bams
        # # # (these two commands are joined because I want to point only
        # # # to the target of the second)
        cmd = """samtools merge \\\n{0} \\\n{1};\n\nsambamba sort -t 8 \\\n{0}\n""".format(
            output_bam, " \\\n".join(bams)
        )
        cmds += [add_cmd(cmd, target=output_sorted_bam, overwrite=overwrite)]

        # # bigWig file
        cmd = bam_to_bigwig(
            output_sorted_bam,
            output_bigwig,
            genome,
            normalization_method=normalization_method,
            cpus=cpus) + "\n"
        cmds += [add_cmd(cmd, target=output_bigwig, overwrite=overwrite)]

        if nucleosome:
            cmd = (
                "sambamba view -f bam \\\n"
                "-t {} \\\n-o {} \\\n-F '(template_length < 100) "
                "and (template_length > -100)' \\\n{}\n"
            ).format(cpus, output_nucleosome_free_reads, output_sorted_bam)
            cmds += [
                add_cmd(cmd, target=output_nucleosome_free_reads, overwrite=overwrite)
            ]
            cmd = (
                "sambamba view -f bam -t {} \\\n-o {} \\\n"
                "-F '((template_length > 180) and (template_length < 247)) "
                "or ((template_length < -180) and (template_length > -247))' \\\n{}\n"
            ).format(cpus, output_nucleosome_reads, output_sorted_bam)
            cmds += [add_cmd(cmd, target=output_nucleosome_reads, overwrite=overwrite)]

            cmd = bam_to_bigwig(
                output_nucleosome_reads,
                output_nucleosome_reads_bigwig, genome) + "\n"
            cmds += [add_cmd(
                cmd,
                target=output_nucleosome_reads_bigwig,
                overwrite=overwrite)]

            cmd = bam_to_bigwig(
                output_nucleosome_free_reads,
                output_nucleosome_free_reads_bigwig, genome) + "\n"
            cmds += [add_cmd(
                cmd,
                target=output_nucleosome_free_reads_bigwig,
                overwrite=overwrite)]

        job = "\n".join(cmds + ["date"]) + "\n"

        with open(job_file, "w") as handle:
            handle.write(job)

        if dry_run:
            continue
        if as_job:
            subprocess.call(
                (
                    "sbatch -p longq --time 4-00:00:00 "
                    "-J merge_signal.{} -o {} -c {} --mem 80000 {}")
                .format(name, log_file, cpus, job_file)
                .split(" ")
            )
        else:
            subprocess.call("sh {}".format(job_file).split(" "))


def bam_to_bigwig(
        input_bam, output_bigwig, genome_assembly,
        normalization_method="RPGC", cpus="8"):
    """
    Convert BAM file to BigWig format using deeptools bamCoverage.

    Parameters
    ----------
    input_bam : str
        Input BAM file.
    output_bigwig : str
        Output BigWig file.
    genome_assembly : str
        Genome assembly of the BAM file.
        One of 'hg19', 'hg38', 'mm10', 'mm9'.
        Otherwise will assume size of human genome.
    normalization_method : {str}, optional
        Method to normalize tracks regarding sequenced depth.
        One of the methods in https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html#Read%20coverage%20normalization%20options

        Default is "RPGC".

    Returns
    -------
    str
        Command to execute.
    """
    from collections import defaultdict

    if genome_assembly not in ['hg19', 'hg38', 'mm10', 'mm9']:
        print(
            "Genome assembly is not known."
            "Using size of human genome. Beware.")

    genome_size = defaultdict(lambda: 3300000000)
    for gen in ['mm9', 'mm10']:
        genome_size[gen] = 2800000000

    cmd = "bamCoverage \\\n--bam {bam_file} \\\n-o {bigwig} \\\n"
    cmd += "-p {cpus} --binSize 10  --normalizeUsing {norm} \\\n"
    cmd += "--effectiveGenomeSize {genome_size} --extendReads 175"""
    cmd = cmd.format(
        bam_file=input_bam,
        bigwig=output_bigwig,
        norm=normalization_method,
        genome_size=genome_size[genome_assembly],
        cpus=cpus)
    return cmd


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
