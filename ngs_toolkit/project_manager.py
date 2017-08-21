#!/usr/bin/env python

"""
A project creator. https://github.com/afrendeiro/tookit
"""

import argparse
import os
import sys
import textwrap


def parse_arguments():
    """
    Argument Parsing.
    """
    description = "%(prog)s - A project creator."
    epilog = "https://github.com/afrendeiro/toolkit"

    parser = argparse.ArgumentParser(
        description=description, epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # parser.add_argument("-V", "--version", action="version",
    #               version="%(prog)s {v}".format(v=__version__))

    parser.add_argument(
        dest="project_name",
        help="Project name.")
    parser.add_argument(
        "-d",
        "--dry-run",
        action="store_true",
        help="Don't actually do anything.")

    # To enable the loop to pass args directly on to the pipelines...
    args = parser.parse_args()

    return args


def create_project(
        project_name, overwrite=False,
        username="arendeiro", email="arendeiro@cemm.oeaw.ac.at",
        url="http://biomedical-sequencing.at/bocklab/arendeiro/{project_name}"):
    """
    Main function: Create project.
    """
    # Easier change later, especially likely for library --> protocol.
    project_dir = os.path.join(os.path.curdir, project_name)

    if os.path.exists(project_dir):
        if not overwrite:
            return

    metadata_dir = os.path.join(project_dir, "metadata")
    project_config = os.path.join(metadata_dir, "project_config.yaml")
    annotation_table = os.path.join(metadata_dir, "annotation.csv")
    merge_table = os.path.join(metadata_dir, "merge_table.csv")
    comparison_table = os.path.join(metadata_dir, "comparison_table.csv")
    src_dir = os.path.join(project_dir, "src")

    # make dirs
    os.makedirs(project_dir)
    os.makedirs(metadata_dir)
    os.makedirs(src_dir)

    if "{project_name}" in url:
        url = url.format(project_name=project_name)

    project_config_template = """    metadata:
        output_dir: /scratch/lab_bock/shared/projects/{project_name}
        results_subdir: data
        submission_subdir: submission
        pipelines_dir: /home/arendeiro/workspace/pipelines
        sample_annotation: /scratch/lab_bock/shared/projects/{project_name}/metadata/annotation.csv
        merge_table: /scratch/lab_bock/shared/projects/{project_name}/metadata/merge_table.csv
    data_sources:
        local: "/scratch/users/arendeiro/data/external/atac-seq/{{sample_name}}.bam"
        bsf: /scratch/lab_bsf/samples/{{flowcell}}/{{flowcell}}_{{lane}}_samples/{{flowcell}}_{{lane}}#{{BSF_name}}.bam
    genomes:
        human: hg19
        mouse: mm10
    transcriptomes:
        human: hg19_cdna
        mouse: mm10_cdna
    pipeline_config:
        atacseq: null
    compute:
        submission_template: templates/slurm_template.sub
        submission_command: sbatch
    trackhubs:
    trackhub_dir: /data/groups/lab_bock/public_html/arendeiro/{project_name}/
    url: {url}
    username: {username}
    email: {email}""".format(project_name=project_name, username=username, email=email, url=url)

    merge_table_template = ",".join([
        "sample_name", "flowcell", "lane", "BSF_name", "data_source"])
    annotation_table_template = ",".join([
        "sample_name", "toggle", "pass_qc", "library",
        "cell_line", "cell_type", "condition",
        "experimental_batch", "experiment_name", "replicate",
        "organism", "flowcell", "lane", "BSF_name", "data_source"])
    comparison_table_template = ",".join([
        "data_type", "comparison_name", "comparison_side",
        "sample_name", "sample_group", "toggle"])

    # write config and tables
    with open(project_config, "w") as handle:
        handle.write(textwrap.dedent(project_config_template))
    with open(merge_table, "w") as handle:
        handle.write(merge_table_template)
    with open(annotation_table, "w") as handle:
        handle.write(annotation_table_template)
    with open(comparison_table, "w") as handle:
        handle.write(comparison_table_template)

    # Initialize git repository
    return os.system("cd {}; git init".format(project_name))


def main():
    """
    Program's main entry point.
    """
    # Parse command-line arguments.
    args = parse_arguments()

    # Create project.
    create_project(args.project_name)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
