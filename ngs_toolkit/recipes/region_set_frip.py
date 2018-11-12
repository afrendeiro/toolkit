#!/usr/bin/env python

"""
This is the "region_set_frip" recipe for ngs_toolkit.
"""


from argparse import ArgumentParser
import os
import sys

from peppy import Project
from ngs_toolkit.atacseq import ATACSeqAnalysis
from ngs_toolkit.chipseq import ChIPSeqAnalysis


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
        "-d", "--data-type",
        dest="data_type",
        default=None,
        help="Data types to perform analysis on. Will be done separately for each.",
        type=str)
    parser.add_argument(
        "-n", "--analysis-name",
        dest="name",
        default=None,
        help="Name of analysis. Will be the prefix of output_files. "
             "By default it will be the name of the Project given in the YAML configuration.",
        type=str)
    parser.add_argument(
        "-r", "--region-set",
        dest="region_set",
        default=None,
        help="BED file with region set derived from several samples or Oracle region set. "
        "If unset, will try to get the `sites` attribute of an existing analysis object "
        "if existing, otherwise will create a region set from the peaks of all samples.",
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
        "-o", "--results-output",
        default="results",
        dest="output_dir",
        help="Directory for analysis output files. "
        "Default is 'results' under the project roort directory.",
        type=str)
    parser.add_argument(
        "-s", "--strict",
        action="store_true",
        dest="strict",
        help="Whether to throw an error in case files cannot be created or not.")
    return parser


def main():
    parser = ArgumentParser(
        prog="region_set_frip",
        description="Region set FRiP recipe."
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

    # ANALYSIS
    if args.data_type is None:
        print("Type of analysis not specified. Will run independent analysis for all types of data in the sample annotation sheet.")
        data_types = sorted(list(set([s.library for s in prj.samples])))
        print("Sample data types: '{}'.".format(",".join(data_types)))
    else:
        print("Type of analysis specified. Will run only analysis for samples of type '{}'.".format(args.data_type))
        data_types = [args.data_type]
        print("Sample data types: '{}'.".format(",".join(data_types)))
    if args.name is None:
        print("Analysis name not specified, will use name in project configuration file: '{}'.".format(prj.project_name))
        args.name = prj.project_name

    for data_type in data_types:
        print("Starting analysis for samples of type: '{}'.".format(data_type))
        samples = [s for s in prj.samples if (s.library == data_type)]
        if len(samples) > 0:
            print(
                "Samples under consideration: '{}'. ".format(",".join([s.name for s in samples])) +
                "Total of {} samples.".format(len([s.name for s in samples])))
        else:
            raise ValueError("There were no valid samples for this analysis type!")

        if data_type in ["ATAC-seq"]:
            print("Initializing ATAC-seq analysis")
            analysis = ATACSeqAnalysis(
                name=args.name + "_atacseq", prj=prj,
                samples=samples, results_dir=args.output_dir)
        elif data_type in ["ChIP-seq", "ChIPmentation"]:
            print("Initializing ChIP-seq analysis")
            analysis = ChIPSeqAnalysis(
                name=args.name + "_chipseq", prj=prj,
                samples=samples, results_dir=args.output_dir)

        region_set_frip = None
        if args.region_set is not None:
            print("Using given region set in BED format: '{}'".format(args.region_set))
            # Use given region set if passed
            region_set_frip = calculate_region_set_frip(
                analysis, region_set=args.region_set, samples=samples,
                as_job=args.as_job, permissive=not args.strict)
        else:
            # try to load an analysis object that contains a `sites` attribute
            print("Trying to load an existing analysis object from pickle file: '{}'".format(analysis.pickle))
            try:
                analysis = analysis.from_pickle()
            except IOError:
                print("Couldn't load an existing analysis object")
                continue

            # if it has, use it
            if hasattr(analysis, "sites"):
                print("Using existing regions set from analysis. Contains {} regions.".format(len(analysis.sites)))
                region_set_frip = calculate_region_set_frip(
                    analysis, region_set=analysis.sites, samples=samples,
                    as_job=args.as_job, permissive=not args.strict)
            else:
                print("Existing analysis does not have an existing region set.")

        if region_set_frip is None:
            print("Generating a new region set for this analysis.")
            # if none worked, make one
            analysis.get_consensus_sites()

            # and use it
            print("Using generated new region set for this analysis.")
            region_set_frip = calculate_region_set_frip(
                analysis, region_set=analysis.sites, samples=samples,
                as_job=args.as_job, permissive=not args.strict)


def calculate_region_set_frip(
        analysis,
        region_set=None, samples=None,
        as_job=False, cpus=8, permissive=True):
    """
    """
    if region_set is None:
        region_set = analysis.sites

    if samples is None:
        samples = analysis.samples

    for sample in samples:
        inside_reads = os.path.join(sample.paths.sample_root, "region_set_frip.inside_reads.txt")
        all_reads = os.path.join(sample.paths.sample_root, "region_set_frip.all_reads.txt")

        log_file = os.path.join(sample.paths.sample_root, "region_set_frip.log")
        job_file = os.path.join(sample.paths.sample_root, "region_set_frip.sh")
        sample_stats = os.path.join(sample.paths.sample_root, "stats.tsv")

        job = "\n".join([
            "#!/bin/sh",
            "date",
            """samtools view -@ {} -c -L {} {} > {}""".format(cpus, region_set, sample.filtered, inside_reads),
            """samtools view -@ {} -c {} > {}""".format(cpus, sample.filtered, all_reads),
            'calc(){ awk "BEGIN { print "$*" }"; }',
            "IN=`cat {}`".format(inside_reads),
            "ALL=`cat {}`".format(all_reads),
            "FRIP=`calc $IN/$ALL`",
            'echo -e "region_set_frip\\t$FRIP\\t." >> {}'.format(sample_stats),
            "date"])

        try:
            with open(job_file, "w") as handle:
                handle.write(job)
        except IOError as e:
            if permissive:
                print("Couldn't write job for sample '{}'.".format(sample.name))
                continue
            else:
                raise e

        if as_job:
            os.system("sbatch -p shortq -J region_set_frip.{} -o {} -c {} --mem 8000 {}".format(sample.name, log_file, cpus, job_file))
        else:
            os.system("sh {}".format(job_file))

    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
