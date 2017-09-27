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
    parser.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help="Don't overwrite any existing directory or file.")

    # To enable the loop to pass args directly on to the pipelines...
    args = parser.parse_args()

    return args


def create_project(
        project_name, overwrite=False,
        username="arendeiro", email="{username}@cemm.oeaw.ac.at",
        url="http://biomedical-sequencing.at/bocklab/{username}/{project_name}"):
    """
    Main function: Create project.
    """
    project_dir = os.path.join(os.path.curdir, project_name)

    if os.path.exists(project_dir):
        if not overwrite:
            print("Detected existing project directory, skipping.")
            return 1

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
        url = url.format(username=username, project_name=project_name)

    if "{username}" in email:
        email = email.format(username=username)


    project_config_template = """    project_name: {project_name}
    project_description: {project_name}
    username: {username}
    email: {email}
    metadata:
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
        url: {url}""".format(
            project_name=project_name, username=username, email=email, url=url)

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
        handle.write(textwrap.dedent(project_config_template + "\n"))
    with open(merge_table, "w") as handle:
        handle.write(merge_table_template)
    with open(annotation_table, "w") as handle:
        handle.write(annotation_table_template)
    with open(comparison_table, "w") as handle:
        handle.write(comparison_table_template)

    # Initialize git repository
    return os.system("cd {}; git init".format(project_name))


def create_requirements_file(
        project_name,
        libraries=[
            "numpy", "scipy", "pandas",
            "matplotlib", "seaborn",
            "pysam", "pybedtools",
            "scikit-learn", "statsmodels", "patsy",
            "looper", "pypiper"],
        overwrite=False):
    """
    Create a requirements.txt file with pip requirements.
    """
    project_dir = os.path.join(os.path.curdir, project_name)
    requirements_file = os.path.join(project_dir, "requirements.txt")

    if os.path.exists(requirements_file):
        if not overwrite:
            print("Detected existing, skipping.")
            return

    requirements_filecontent = "\n".join(libraries)

    # write requirements file
    with open(requirements_file, "w") as handle:
        handle.write(textwrap.dedent(requirements_filecontent) + "\n")


def create_makefile(
        project_name, overwrite=False):
    """
    Create a Makefile to manage the project execution.
    """
    project_dir = os.path.join(os.path.curdir, project_name)
    makefile = os.path.join(project_dir, "Makefile")
    src_dir = os.path.join(project_dir, "src")
    log_dir = os.path.join(project_dir, "log")
    metadata_dir = os.path.join(project_dir, "metadata")
    project_config = os.path.join(metadata_dir, "project_config.yaml")

    if os.path.exists(makefile):
        if not overwrite:
            print("Detected existing, skipping.")
            return

    makefile_content = """    .DEFAULT_GOAL := analysis_job

    requirements:
        pip install -r requirements.txt

    # process project's data
    # with looper/pypiper/pipelines:
    # see https://github.com/epigen/looper
    # see https://github.com/epigen/pypiper
    # see https://github.com/epigen/open_pipelines
    process:
        looper run {project_config}

    analysis:
        looper summarize {project_config}
        python -u src/analysis.py

    analysis_job:
        sbatch -p shortq -c 12 --mem 80000 -J {project_name}.analysis -o {log_dir}/{project_name}.analysis.log --wrap "python -u src/analysis.py"

    all: requirements process analysis

    .PHONY: requirements process analysis all""".format(
        project_config=project_config,
        project_name=project_name,
        log_dir=log_dir)

    # write Makefile
    with open(makefile, "w") as handle:
        handle.write(textwrap.dedent(makefile_content) + "\n")


def main():
    """
    Program's main entry point.
    """
    # Parse command-line arguments.
    args = parse_arguments()

    # Create project.
    git_ok = create_project(
        project_name=args.project_name,
        overwrite=args.overwrite)
    if git_ok != 0:
        return git_ok

    # Create requirements file.
    create_requirements_file(
        project_name=args.project_name,
        overwrite=args.overwrite)

    # Create Makefile.
    create_makefile(
        project_name=args.project_name,
        overwrite=args.overwrite)

    


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
