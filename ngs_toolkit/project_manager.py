#!/usr/bin/env python

"""
A project creator. https://github.com/afrendeiro/tookit
"""

import argparse
import os
import sys
import textwrap
import ngs_toolkit


_LOGGER = ngs_toolkit._LOGGER


def parse_arguments():
    """
    Argument Parsing.
    """
    description = "%(prog)s - A project manager."
    epilog = "https://github.com/afrendeiro/toolkit"

    parser = argparse.ArgumentParser(
        description=description, epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-V", "--version", action="version",
                        version=ngs_toolkit.__version__)

    subparsers = parser.add_subparsers(dest="command")

    # Create command
    create_subparser = subparsers.add_parser(
        "create", description="Create project.",
        help="Create project.")
    create_subparser.add_argument(
        dest="project_name",
        help="Project name.")
    create_subparser.add_argument(
        '-r', '--root-dir',
        default=os.path.curdir,
        dest="root_dir",
        help="Root directory to create projects.")
    create_subparser.add_argument(
        "-d",
        "--dry-run",
        action="store_true",
        help="Don't actually do anything.")
    create_subparser.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help="Don't overwrite any existing directory or file.")

    # Recipe command
    recipe_subparser = subparsers.add_parser(
        "recipe", description="Run recipe.",
        help="Run ngs_toolkit recipe for a given project.")
    recipe_subparser.add_argument(
        dest="recipe_name",
        help="Recipe name.")
    recipe_subparser.add_argument(
        dest="project_config",
        help="Project configuration file.")

    for p in [create_subparser, recipe_subparser]:
        p.add_argument("-V", "--version", action="version",
                       version=ngs_toolkit.__version__)

    args = parser.parse_args()
    args.root_dir = os.path.abspath(args.root_dir)

    return args


def create_project(
        project_name, root_dir, overwrite=False,
        username="arendeiro", email="{username}@cemm.oeaw.ac.at",
        url="http://biomedical-sequencing.at/bocklab/{username}/{project_name}"):
    """
    Main function: Create project.
    """
    project_dir = os.path.join(root_dir, project_name)

    if os.path.exists(project_dir):
        if not overwrite:
            _LOGGER.error("Detected existing project directory, skipping.")
            return 1

    metadata_dir = os.path.join(project_dir, "metadata")
    project_config = os.path.join(metadata_dir, "project_config.yaml")
    annotation_table = os.path.join(metadata_dir, "annotation.csv")
    merge_table = os.path.join(metadata_dir, "merge_table.csv")
    comparison_table = os.path.join(metadata_dir, "comparison_table.csv")
    src_dir = os.path.join(project_dir, "src")

    # make dirs
    for d in [project_dir, metadata_dir, src_dir]:
        os.makedirs(d, exist_ok=True)

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
        pipeline_interfaces: /home/{username}/workspace/open_pipelines/pipeline_interface.yaml
        sample_annotation: /scratch/lab_bock/shared/projects/{project_name}/metadata/annotation.csv
        sample_subannotation: /scratch/lab_bock/shared/projects/{project_name}/metadata/merge_table.csv
        comparison_table: /scratch/lab_bock/shared/projects/{project_name}/metadata/comcomparison_table.csv
    sample_attributes:
        - sample_name
    group_attributes:
        - sample_name
    data_sources:
        local: "/scratch/users/{username}/data/external/atac-seq/{{sample_name}}.bam"
        bsf: /scratch/lab_bsf/samples/{{flowcell}}/{{flowcell}}_{{lane}}_samples/{{flowcell}}_{{lane}}#{{BSF_name}}.bam
    implied_columns:
        organism:
            "Homo sapiens":
                genome: hg19
                transcriptome: hg19_cdna
            "human":
                genome: hg19
                transcriptome: hg19_cdna
            "Mus musculus":
                genome: mm10
                transcriptome: mm10_cdna
            "mouse":
                genome: mm10
                transcriptome: mm10_cdna
    pipeline_config:
        atacseq: null
    compute:
        submission_template: slurm_template.sub
        submission_command: sbatch
    trackhubs:
        trackhub_dir: /data/groups/lab_bock/public_html/{username}/{project_name}/
        url: {url}""".format(
        project_name=project_name, username=username, email=email, url=url)

    merge_table_template = ",".join([
        "sample_name", "flowcell", "lane", "BSF_name", "data_source"])
    annotation_table_template = ",".join([
        "sample_name", "toggle", "pass_qc", "protocol", "library",
        "cell_line", "cell_type", "condition",
        "experimental_batch", "experiment_name", "replicate",
        "organism", "flowcell", "lane", "BSF_name", "data_source"])
    comparison_table_template = ",".join([
        "comparison_type", "data_type", "comparison_name", "comparison_side",
        "sample_name", "sample_group", "comparison_genome", "toggle"])

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
    os.chdir(project_dir)
    return os.system("git init")


def create_requirements_file(
        project_name,
        project_dir,
        libraries=["numpy>=1.14.5",
                   "scipy>=1.0.1",
                   "pandas>=0.23.3",
                   "matplotlib>=2.1.2",
                   "seaborn>=0.9.0",
                   "pysam>=0.14.1",
                   "pybedtools>=0.7.10",
                   "scikit-learn>=0.20.0",
                   "statsmodels>=0.9.0",
                   "patsy>=0.5.0",
                   "rpy2==2.8.6",
                   "peppy==v0.9.2",
                   "piper==0.6.0",
                   "ngs-toolkit>={__version__}"],
        overwrite=False):
    """
    Create a requirements.txt file with pip requirements.
    """
    libraries = [l.format(__version__=ngs_toolkit.__version__)
                 if "__version__" in l else l for l in libraries]

    requirements_file = os.path.join(project_dir, "requirements.txt")

    if os.path.exists(requirements_file):
        if not overwrite:
            _LOGGER.warn("'requirements.txt' file already existing, skipping.")
            return

    requirements_filecontent = "\n".join(libraries)

    # write requirements file
    with open(requirements_file, "w") as handle:
        handle.write(textwrap.dedent(requirements_filecontent) + "\n")


def create_makefile(
        project_name, project_dir, overwrite=False):
    """
    Create a Makefile to manage the project execution.
    """
    makefile = os.path.join(project_dir, "Makefile")
    src_dir = "src"
    log_dir = "log"
    metadata_dir = "metadata"
    project_config = os.path.join(metadata_dir, "project_config.yaml")

    if os.path.exists(makefile):
        if not overwrite:
            print("Detected existing, skipping.")
            return

    makefile_content = """    .DEFAULT_GOAL := analysis_job

    requirements:
        pip install -r requirements.txt

    process:
        looper run {project_config}

    summarize:
        looper summarize {project_config}

    mklog:
        mkdir -p log

    analysis: summarize
        ngs_analysis.py {project_config}

    analysis_job: summarize mklog
        sbatch -p longq --time 8-00:00:00 -c 12 --mem 80000 \\
        -J {project_name}.analysis \\
        -o {log_dir}/$(shell date +"%Y%m%d-%H%M%S").{project_name}.analysis.log \\
        --x11 --wrap "ngs_analysis {project_config}"

    all: requirements process analysis

    .PHONY: requirements process summarize mklog analysis all""".format(
        project_config=project_config,
        project_name=project_name,
        log_dir=log_dir).replace("    ", "\t")

    # write Makefile
    with open(makefile, "w") as handle:
        handle.write(textwrap.dedent(makefile_content) + "\n")


def run_recipe(recipe_name, project_config):
    return os.popen("python -m ngs_toolkit.recipes.{} {}".format(recipe_name, project_config))


def main():
    """
    Program's main entry point.
    """
    # Parse command-line arguments.
    _LOGGER.debug("Parsing command-line arguments")
    args = parse_arguments()

    if args.command == "create":
        _LOGGER.info("Creating project '{}' in '{}'.".format(
            args.project_name, args.root_dir))
        # Create project.
        git_ok = create_project(
            project_name=args.project_name,
            root_dir=args.root_dir,
            overwrite=args.overwrite)
        if git_ok != 0:
            _LOGGER.error("Initialization of project failed.")
            return git_ok

        # Create requirements file.
        _LOGGER.info("Creating requirements file for project '{}'.".format(
            args.project_name))
        create_requirements_file(
            project_name=args.project_name,
            project_dir=os.path.join(args.root_dir, args.project_name),
            overwrite=args.overwrite)

        # Create Makefile.
        _LOGGER.info("Creating Makefile file for project '{}'.".format(
            args.project_name))
        create_makefile(
            project_name=args.project_name,
            project_dir=os.path.join(args.root_dir, args.project_name),
            overwrite=args.overwrite)

    elif args.command == "recipe":
        _LOGGER.info("Running recipe '{}'.".format(args.recipe_name))
        run_recipe(
            recipe_name=args.recipe_name,
            project_config=args.project_config)

    _LOGGER.info("Completed.")


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
