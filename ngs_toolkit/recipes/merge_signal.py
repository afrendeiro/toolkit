#!/usr/bin/env python

"""
This is the "merge_signal" recipe for ngs_toolkit.
"""


from argparse import ArgumentParser
import os
import sys

import pandas as pd

from peppy import Project


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
        "-n", "--normalize",
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
    return parser


def main():
    parser = ArgumentParser(
        prog="merge_signal",
        description="Merge signal recipe."
    )
    parser = add_args(parser)
    args = parser.parse_args()
    # args = parser.parse_args('-t ATAC-seq metadata/project_config.yaml'.split(" "))

    # Start project
    print("Starting peppy project with project configuration file: '{}'".format(args.config_file))
    prj = Project(args.config_file)
    print("Changing directory to project root directory: '{}'.".format(prj.metadata.output_dir))
    os.chdir(prj.metadata.output_dir)
    if args.pass_qc:
        print("Filtering samples out which didn't pass QC as specified in sample annotation in column 'pass_qc'")
        prj._samples = [s for s in prj._samples if s.pass_qc not in ['0', 0, 'False', False]]
    print("Setting location of sample files dependent on sample types.")
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
            "Total of {} samples.".format(len([s.name for s in prj.samples])))
    else:
        raise ValueError("There were no valid samples after filtering for quality!")

    # Get only samples with signal
    key = "protocol" if 'protocol' in prj.sheet.columns else 'library'
    sheet = prj.sheet[prj.sheet[key].isin(["ATAC-seq", "ChIP-seq", "ChIPmentation"])]
    print(
        "Selecting samples with appropriate data type." +
        "Samples under consideration: '{}'. ".format(",".join(sheet['sample_name'].tolist())) +
        "Total of {} samples.".format(sheet.shape[0]))

    # Get default attributes if not set
    if args.attributes is None:
        if "group_attributes" in prj:
            args.attributes = prj.group_attributes
        else:
            print("Sample attributes to group by were not set and none could be found in project configuration file!")
            print("Aborting!")
            return
    else:
        args.attributes = args.attributes.split(",")

    print("Using the following attributes to merge samples: '{}', resulting in a total of {} groups.".format(
        ", ".join(args.attributes), len(sheet.groupby(args.attributes).groups.items())))

    print("Starting to merge sample signals.")
    merge_signal(
        sheet, prj.samples, args.attributes, output_dir=args.output_dir,
        normalize=args.normalize, nucleosome=args.nucleosome, overwrite=args.overwrite, cpus=args.cpus, as_job=args.as_job)


def merge_signal(
        sheet, samples, attributes,
        output_dir="merged",
        normalize=True,
        nucleosome=False,
        overwrite=False,
        cpus=8, as_job=False):
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
        name = "_".join([a for a in attrs if not pd.isnull(a)])
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

        genome_sizes = "/data/groups/lab_bock/shared/resources/genomes/{g}/{g}.chromSizes".format(g=genome)
        transient_file = os.path.abspath(re.sub("\.bigWig", "", output_bigwig))

        # Prepare commands
        cmds = ["#!/bin/sh", "date"]

        ## merge bams
        cmd = """samtools merge {0} {1}; sambamba sort -t 8 {0}""".format(output_bam, " ".join(bams))
        cmds += [add_cmd(cmd, target=output_sorted_bam, overwrite=overwrite)]

        ## bigWig file
        cmd = "; ".join([
            "bedtools bamtobed -i {0} |".format(output_sorted_bam) +
            " bedtools slop -i stdin -g {0} -s -l 0 -r 130 |".format(genome_sizes) +
            " fix_bedfile_genome_boundaries.py {0} |".format(genome) +
            " genomeCoverageBed -bg -g {0} -i stdin > {1}.cov".format(genome_sizes, transient_file),
            "awk 'NR==FNR{{sum+= $4; next}}{{ $4 = ($4 / sum) * 1000000; print}}' {0}.cov {0}.cov > {0}.normalized.cov".format(transient_file) if normalize else "",
            "bedGraphToBigWig {0}{1}.cov {2} {3}".format(transient_file, ".normalized" if normalize else "", genome_sizes, output_bigwig),
            "if [[ -s {0}.cov ]]; then rm {0}.cov; fi".format(transient_file),
            "if [[ -s {0}.normalized.cov ]]; then rm {0}.normalized.cov; fi".format(transient_file)])
        cmds += [add_cmd(cmd, target=output_bigwig, overwrite=overwrite)]

        if nucleosome:
            cmd = ("""sambamba view -f bam -t {} -o {} -F "(template_length < 100) and (template_length > -100)" {}"""
                .format(cpus, output_nucleosome_free_reads, output_sorted_bam))
            cmds += [add_cmd(cmd, target=output_nucleosome_free_reads, overwrite=overwrite)]
            cmd = ("""sambamba view -f bam -t {} -o {} -F "((template_length > 180) and (template_length < 247)) or ((template_length < -180) and (template_length > -247))" {}"""
                .format(cpus, output_nucleosome_reads, output_sorted_bam))
            cmds += [add_cmd(cmd, target=output_nucleosome_reads, overwrite=overwrite)]

        job = "\n".join(cmds + ['date'])

        with open(job_file, "w") as handle:
            handle.write(job)

        if as_job:
            os.system("sbatch -p longq --time 4-00:00:00 -J merge_signal.{} -o {} -c {} --mem 80000 {}".format(name, log_file, cpus, job_file))
        else:
            os.system("sh {}".format(job_file))


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
