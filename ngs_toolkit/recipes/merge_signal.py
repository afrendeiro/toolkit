#!/usr/bin/env python

"""
This is the "merge_signal" recipe for ngs_toolkit.
"""


from argparse import ArgumentParser
import os
import sys

import pandas as pd

from peppy import Project
from ngs_toolkit import _LOGGER, _CONFIG, __version__


class Unbuffered(object):
    def __init__(self, stream):
        self.stream = stream

    def write(self, data):
        self.stream.write(data)
        self.stream.flush()

    def writelines(self, datas):
        self.stream.writelines(datas)
        self.stream.flush()

    def __getattr__(self, attr):
        return getattr(self.stream, attr)


sys.stdout = Unbuffered(sys.stdout)


def add_args(parser):
    """
    Global options for analysis.
    """
    parser.add_argument(
        dest="config_file",
        help="YAML project configuration file.",
        type=str)
    parser.add_argument(
        "-a", "--attributes",
        dest="attributes",
        default=None,
        help="Attributes to merge samples by. A comma-delimited string with no spaces. "
             "By default will use values in the project config `group_attributes`.",
        type=str)
    parser.add_argument(
        "-q", "--pass-qc",
        action="store_true",
        dest="pass_qc",
        help="Whether only samples with a 'pass_qc' value of '1' "
        "in the annotation sheet should be used.")
    parser.add_argument(
        "-j", "--as-jobs",
        action="store_true",
        dest="as_job",
        help="Whether jobs should be created for each sample, or "
        "it should run in serial mode.")
    parser.add_argument(
        "--cpus",
        dest="cpus",
        default=8,
        help="CPUs/Threads to use per job if `--as-jobs` is on.")
    parser.add_argument(
        "--normalize",
        action="store_true",
        dest="normalize",
        help="Whether tracks should be normalized to total sequenced "
        "depth.")
    parser.add_argument(
        "--nucleosome",
        action="store_true",
        dest="nucleosome",
        help="Whether to produce nucleosome/nucleosome-free signal files.")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        dest="overwrite",
        help="Whether to overwrite existing files.")
    parser.add_argument(
        "-o", "--output-dir",
        default="merged",
        dest="output_dir",
        help="Directory for output files. "
        "Default is 'merged' under the project root directory.",
        type=str)
    parser.add_argument(
        "-d", "--dry-run",
        action="store_true",
        dest="dry_run",
        help="Whether to do everything except running commands.")
    return parser


def main():
    desc = """
    Merge signal recipe.
    This recipe will merge signal from various ATAC-seq or ChIP-seq samples \
    given a set of attributes to group samples by.

    It produces merged BAM and bigWig files for all signal in the samples but \
    is also capable of producing this for nucleosomal/nucleosomal free signal \
    based on fragment length distribution if data is paired-end sequenced. \
    This signal may optionally be normalized for each group. \
    It is also capable of parallelizing work in jobs if a SLURM cluster is available."""

    parser = ArgumentParser(prog="merge_signal", description=desc)
    parser = add_args(parser)
    args = parser.parse_args()
    # args = parser.parse_args('-t ATAC-seq metadata/project_config.yaml'.split(" "))

    _LOGGER.info("This is the 'merge_signal' recipe from ngs_toolkit, version: {}".format(__version__))
    # Start project
    _LOGGER.debug("Starting peppy project with project configuration file: '{}'".format(args.config_file))
    prj = Project(args.config_file)
    _LOGGER.debug("Changing directory to project root directory: '{}'.".format(prj.metadata.output_dir))
    os.chdir(prj.metadata.output_dir)
    if args.pass_qc:
        _LOGGER.info("Filtering samples out which didn't pass QC as specified in sample \
                      annotation in column 'pass_qc'")
        prj._samples = [s for s in prj._samples if s.pass_qc not in ['0', 0, 'False', False]]
    _LOGGER.debug("Setting location of sample files dependent on sample types.")
    for sample in prj.samples:
        if hasattr(sample, "protocol"):
            sample.library = sample.protocol

        if sample.library in ["ATAC-seq", "ChIP-seq", "ChIPmentation"]:
            sample.mapped = os.path.join(
                sample.paths.sample_root,
                "mapped", sample.name + ".trimmed.bowtie2.bam")
            sample.filtered = os.path.join(
                sample.paths.sample_root,
                "mapped", sample.name + ".trimmed.bowtie2.filtered.bam")
            sample.peaks = os.path.join(
                sample.paths.sample_root,
                "peaks", sample.name + "_peaks.narrowPeak")

    if len(prj.samples) > 0:
        print(
            "Samples under consideration: '{}'. ".format(",".join([s.name for s in prj.samples])) +
            "\nTotal of {} samples.".format(len([s.name for s in prj.samples])))
    else:
        raise ValueError("There were no valid samples after filtering for quality!")

    # Get only samples with signal
    key = "protocol" if 'protocol' in prj.sheet.columns else 'library'
    sheet = prj.sheet[prj.sheet[key].isin(["ATAC-seq", "ChIP-seq", "ChIPmentation"])]
    sheet = sheet.loc[sheet['sample_name'].isin([s.name for s in prj.samples])]
    _LOGGER.info(
        "Selecting samples with appropriate data type." +
        "\nSamples under consideration: '{}'. ".format(",".join(sheet['sample_name'].tolist())) +
        "\nTotal of {} samples.".format(sheet.shape[0]))

    # Get default attributes if not set
    if args.attributes is None:
        if "group_attributes" in prj:
            args.attributes = prj.group_attributes
        else:
            _LOGGER.error("Sample attributes to group by were not set and none could be found in project \
                          configuration file!" + "\nAborting!")
            return 1
    else:
        if "," in args.attributes:
            args.attributes = args.attributes.split(",")
        else:
            args.attributes = [args.attributes]

    _LOGGER.info("Using the following attributes to merge samples: '{}', resulting in a total of {} groups."
                 .format(", ".join(args.attributes), len(sheet.groupby(args.attributes).groups.items())))

    merge_signal(
        sheet, prj.samples, args.attributes, output_dir=args.output_dir,
        normalize=args.normalize, nucleosome=args.nucleosome,
        overwrite=args.overwrite, cpus=args.cpus, as_job=args.as_job, dry_run=args.dry_run)


def merge_signal(
        sheet, samples, attributes,
        output_dir="merged",
        normalize=True,
        nucleosome=False,
        overwrite=False,
        cpus=8, as_job=False, dry_run=False):
    """
    """
    import re

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
        if isinstance(attrs, list):
            name = "_".join([a for a in attrs if not pd.isnull(a)])
        else:
            name = attrs
        name = re.sub(r"_{2}", "_", name)
        name = re.sub(r"_$", "", name)
        sample_subset = [s for s in samples if s.name in sheet.loc[index, "sample_name"].tolist()]
        bams = [s.filtered for s in sample_subset]

        genomes = list(set(s.genome for s in sample_subset))
        if len(genomes) != 1:
            raise ValueError("Samples under same attribute group have more than one genome assembly: '{}'.".format("', '".join(genomes)))
        else:
            genome = genomes[0]

        job_file = os.path.join(output_dir, name + ".merge_signal.sh")
        log_file = os.path.join(output_dir, name + ".merge_signal.log")
        output_bam = os.path.join(output_dir, name + ".merged.bam")
        output_sorted_bam = os.path.join(output_dir, name + ".merged.sorted.bam")
        output_bigwig = os.path.join(output_dir, name + ".bigWig")
        output_nucleosome_free_reads = os.path.join(output_dir, name + ".nucleosome_free_reads.bam")
        output_nucleosome_reads = os.path.join(output_dir, name + ".nucleosome_free_reads.bam")

        # Get region databases from config
        _LOGGER.debug("Getting genome assembly chromosome sizes for genome '{}' from configuration.".format(genome))

        msg = "Genome assembly chromosome sizes database values in configuration could not be found or understood. "
        msg += "Please add a path to this file to this section 'resources:genomes:chrom_sizes:{}'. ".format(genome)
        msg += "For an example, see https://github.com/afrendeiro/toolkit/tree/master/ngs_toolkit/config/example.yaml"
        try:
            chrom_sizes = _CONFIG['resources']['genomes']['chrom_sizes'][genome]
        except KeyError:
            _LOGGER.error(msg)
            return
        transient_file = os.path.abspath(re.sub(r"\.bigWig", "", output_bigwig))

        # Prepare commands
        cmds = ["#!/bin/sh", "date"]

        # # merge bams
        cmd = """samtools merge {0} {1}; sambamba sort -t 8 {0}""".format(output_bam, " ".join(bams))
        cmds += [add_cmd(cmd, target=output_sorted_bam, overwrite=overwrite)]

        # # bigWig file
        cmd = "\n".join([
            "bedtools bamtobed -i {0} |".format(output_sorted_bam) +
            " bedtools slop -i stdin -g {0} -s -l 0 -r 130 |".format(chrom_sizes) +
            " fix_bedfile_genome_boundaries.py {0} |".format(genome) +
            " genomeCoverageBed -bg -g {0} -i stdin > {1}.cov".format(chrom_sizes, transient_file),
            "awk 'NR==FNR{{sum+= $4; next}}{{ $4 = ($4 / sum) * 1000000; print}}' {0}.cov {0}.cov > {0}.normalized.cov".format(transient_file) if normalize else "",
            "bedGraphToBigWig {0}{1}.cov {2} {3}".format(transient_file, ".normalized" if normalize else "", chrom_sizes, output_bigwig),
            "if [[ -s {0}.cov ]]\nthen rm {0}.cov\nfi".format(transient_file),
            "if [[ -s {0}.normalized.cov ]]\nthen rm {0}.normalized.cov\nfi".format(transient_file)])
        cmds += [add_cmd(cmd, target=output_bigwig, overwrite=overwrite)]

        if nucleosome:
            cmd = ("""sambamba view -f bam -t {} -o {} -F "(template_length < 100) and (template_length > -100)" {}"""
                   .format(cpus, output_nucleosome_free_reads, output_sorted_bam))
            cmds += [add_cmd(cmd, target=output_nucleosome_free_reads, overwrite=overwrite)]
            cmd = ("""sambamba view -f bam -t {} -o {} -F "((template_length > 180) and (template_length < 247)) or ((template_length < -180) and (template_length > -247))" {}"""
                   .format(cpus, output_nucleosome_reads, output_sorted_bam))
            cmds += [add_cmd(cmd, target=output_nucleosome_reads, overwrite=overwrite)]

        job = "\n".join(cmds + ['date']) + "\n"

        with open(job_file, "w") as handle:
            handle.write(job)

        if dry_run:
            continue
        if as_job:
            os.system("sbatch -p longq --time 4-00:00:00 " +
                      "-J merge_signal.{} -o {} -c {} --mem 80000 {}".format(name, log_file, cpus, job_file))
        else:
            os.system("sh {}".format(job_file))


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
