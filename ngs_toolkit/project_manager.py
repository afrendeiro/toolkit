#!/usr/bin/env python

"""
The project manager from ngs_toolkit.

Create and run NGS-based projects with the PEP specification and run analysis using 'ngs_toolkit'.
"""

import argparse
import os
import sys
import textwrap
import ngs_toolkit
from ngs_toolkit import _CONFIG, _LOGGER


def parse_arguments(cli_string=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(
        dest="command"
    )  # , required=True <- not supported in Python < 3.7

    # Create command
    create_subparser = subparsers.add_parser(
        "create", description="Create project.", help="Create project."
    )
    create_subparser.add_argument(dest="project_name", help="Project name.")
    # # parse default assemblies from config
    default = ",".join(
        [
            ":".join([k, v])
            for x in _CONFIG["preferences"]["default_genome_assemblies"]
            for k, v in x.items()
        ]
    )
    # default = "human:hg19,mouse:mm10"
    create_subparser.add_argument(
        "-g",
        "--genome-assembly",
        default=default,
        dest="genome_assemblies",
        help="List of 'organism:assembly' pairs for project. "
        "Comma-separated list of pairs of supported organism/genome assembly. "
        "Defaults to '{}'.".format(default),
    )
    create_subparser.add_argument(
        "-r",
        "--root-dir",
        default=os.path.curdir,
        dest="root_dir",
        help="Root directory to create projects.",
    )
    create_subparser.add_argument(
        "-d", "--dry-run", action="store_true", help="Don't actually do anything."
    )
    create_subparser.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help="Don't overwrite any existing directory or file.",
    )

    # Recipe command
    recipe_subparser = subparsers.add_parser(
        "recipe", description="Run recipe.", help="Run ngs_toolkit recipe for a given project.",
    )
    recipe_subparser.add_argument(dest="recipe_name", help="Recipe name.", nargs="?")
    recipe_subparser.add_argument(
        dest="project_config", help="Project configuration file.", nargs="?"
    )
    recipe_subparser.add_argument(
        "-l",
        "--list",
        dest="list_only",
        action="store_true",
        default=False,
        help="List available recipes and don't do anything else.",
    )

    for p in [create_subparser, recipe_subparser]:
        p.add_argument("-V", "--version", action="version", version=ngs_toolkit.__version__)

    if cli_string is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(cli_string.split(" "))

    if args.command is None:
        parser.print_help(sys.stderr)
        sys.exit(1)
    elif args.command == "create":
        args.root_dir = os.path.abspath(args.root_dir)
    elif args.command == "recipe":
        if not args.list_only and ((args.recipe_name is None) or (args.project_config is None)):
            parser.print_help(sys.stderr)
            sys.exit(1)

    return args


def create_project(
    project_name,
    genome_assemblies,
    overwrite=False,
    root_projects_dir=None,
    username=None,
    email=None,
    url=None,
    git=True,
):
    """
    Main function: Create project.
    """
    import subprocess

    from ngs_toolkit.analysis import Analysis

    # Get defaults from config
    if root_projects_dir is None:
        root_projects_dir = _CONFIG["preferences"]["root_projects_dir"]

    root_projects_dir = Analysis._format_string_with_environment_variables(root_projects_dir)
    project_dir = os.path.join(root_projects_dir, project_name)

    if os.path.exists(project_dir):
        if not overwrite:
            _LOGGER.error("Detected existing project directory, skipping.")
            return 1

    # Get defaults from config
    if username is None:
        username = _CONFIG["username"]
    if username is None:
        username = os.getenv("USER")
    if email is None:
        email = _CONFIG["email"]
    if url is None:
        url = _CONFIG["website_root"]
    if url is not None:
        if "{project_name}" in url:
            url = url.format(project_name=project_name)

    metadata_dir = os.path.join(project_dir, "metadata")
    project_config = os.path.join(metadata_dir, "project_config.yaml")
    annotation_table = os.path.join(metadata_dir, "annotation.csv")
    sample_subannotation = os.path.join(metadata_dir, "sample_subannotation.csv")
    comparison_table = os.path.join(metadata_dir, "comparison_table.csv")
    src_dir = os.path.join(project_dir, "src")

    genome_assemblies = "\n            ".join(
        [
            "- if:{n12}    organism: '{org}'{n12}  then:{n12}    genome: '{gen}'".format(
                org=s, gen=g, n12="\n" + "".join([" "] * 12)
            )
            for s, g in genome_assemblies.items()
        ]
    )

    # make dirs
    for d in [project_dir, metadata_dir, src_dir]:
        if not os.path.exists(d):
            os.makedirs(d)

    project_config_template = """    pep_version: "2.0.0"
    project_name: {project_name}
    description: {project_name}
    username: {username}
    email: {email}
    root_dir: {project_dir}
    results_subdir: data
    submission_subdir: submission
    pipeline_interfaces: /home/{username}/workspace/open_pipelines/pipeline_interface.yaml
    sample_table: {annotation_table}
    subsample_table: {sample_subannotation}
    comparison_table: {comparison_table}
    sample_attributes:
        - sample_name
    group_attributes:
        - sample_name
    sample_modifiers:
        imply:
            {genome_assemblies}
        derive:
            attributes: [data_source]
            sources:
                bsf: /scratch/lab_bsf/samples/{{flowcell}}/{{flowcell}}_{{lane}}_samples/{{flowcell}}_{{lane}}#{{BSF_name}}.bam
                local: /tmp/tmptd4zmpiw/test_project/data/{{sample_name}}.bam
    trackhubs:
        trackhub_dir: {project_dir}/trackhubs
        url: {url}""".format(
        project_name=project_name,
        username=username,
        email=email,
        project_dir=project_dir,
        annotation_table=annotation_table,
        sample_subannotation=sample_subannotation,
        comparison_table=comparison_table,
        genome_assemblies=genome_assemblies,
        url=url,
    )

    merge_table_template = ",".join(["sample_name", "flowcell", "lane", "BSF_name", "data_source"])
    annotation_table_template = ",".join(
        [
            "sample_name",
            "toggle",
            "pass_qc",
            "protocol",
            "library",
            "cell_line",
            "cell_type",
            "condition",
            "experimental_batch",
            "experiment_name",
            "replicate",
            "organism",
            "flowcell",
            "lane",
            "BSF_name",
            "data_source",
        ]
    )
    comparison_table_template = ",".join(
        [
            "comparison_type",
            "data_type",
            "comparison_name",
            "comparison_side",
            "sample_name",
            "sample_group",
            "comparison_genome",
            "toggle",
        ]
    )

    # write config and tables
    with open(project_config, "w", 1) as handle:
        handle.write(textwrap.dedent(project_config_template + "\n"))
    with open(sample_subannotation, "w", 1) as handle:
        handle.write(merge_table_template)
    with open(annotation_table, "w", 1) as handle:
        handle.write(annotation_table_template)
    with open(comparison_table, "w", 1) as handle:
        handle.write(comparison_table_template)

    # Initialize git repository)
    if git:
        p = subprocess.Popen(
            "git init {}".format(project_dir).split(" "),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        p.communicate()
        return p.returncode
    return 0


def create_requirements_file(project_name, project_dir, requirements=None, overwrite=False):
    """
    Create a requirements.txt file with pip requirements.
    """

    def get_current_requirements():
        import requests

        package_name = "ngs-toolkit"
        url = "https://pypi.python.org/pypi/" + str(package_name) + "/json"
        data = requests.get(url).json()
        requirements = [
            x.replace(r" ", "").replace("(", "").replace(")", "")
            for x in data["info"]["requires_dist"]
            if "extra" not in x
        ]
        requirements.append("ngs_toolkit=={}".format(ngs_toolkit.__version__))
        return requirements

    if requirements is None:
        requirements = get_current_requirements()

    requirements_file = os.path.join(project_dir, "requirements.txt")

    if os.path.exists(requirements_file):
        if not overwrite:
            _LOGGER.warning("'requirements.txt' file already existing, skipping.")
            return

    requirements_filecontent = "\n".join(requirements)

    # write requirements file
    with open(requirements_file, "w", 1) as handle:
        handle.write(textwrap.dedent(requirements_filecontent) + "\n")


def create_makefile(project_name, project_dir, overwrite=False):
    """
    Create a Makefile to manage the project execution.
    """
    makefile = os.path.join(project_dir, "Makefile")
    # src_dir = "src"
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
        ngs_analysis {project_config}

    analysis_job: summarize mklog
        sbatch -p longq --time 8-00:00:00 -c 12 --mem 80000 \\
        -J {project_name}.analysis \\
        -o {log_dir}/$(shell date +"%Y%m%d-%H%M%S").{project_name}.analysis.log \\
        --x11 --wrap "ngs_analysis {project_config}"

    all: requirements process analysis

    .PHONY: requirements process summarize mklog analysis all""".format(
        project_config=project_config, project_name=project_name, log_dir=log_dir
    ).replace(
        "    ", "\t"
    )

    # write Makefile
    with open(makefile, "w", 1) as handle:
        handle.write(textwrap.dedent(makefile_content) + "\n")


def run_recipe(recipe_name, project_config):
    import subprocess

    return subprocess.call(
        "{} -m ngs_toolkit.recipes.{} {}".format(sys.executable, recipe_name, project_config).split(
            " "
        )
    )


def main():
    """
    Program"s main entry point.
    """
    # Parse command-line arguments.
    _LOGGER.debug("Parsing command-line arguments")
    args = parse_arguments()

    if args.command == "create":
        _LOGGER.info("Creating project '{}' in '{}'.".format(args.project_name, args.root_dir))

        genome_assemblies = {
            x.split(":")[0]: x.split(":")[1] for x in args.genome_assemblies.split(",")
        }
        # Create project.
        git_ok = create_project(
            project_name=args.project_name,
            genome_assemblies=genome_assemblies,
            overwrite=args.overwrite,
            root_projects_dir=args.root_dir,
        )
        if git_ok != 0:
            _LOGGER.error("Initialization of project failed.")
            return git_ok

        # Create requirements file.
        _LOGGER.info("Creating requirements file for project '{}'.".format(args.project_name))
        create_requirements_file(
            project_name=args.project_name,
            project_dir=os.path.join(args.root_dir, args.project_name),
            overwrite=args.overwrite,
        )

        # Create Makefile.
        _LOGGER.info("Creating Makefile file for project '{}'.".format(args.project_name))
        create_makefile(
            project_name=args.project_name,
            project_dir=os.path.join(args.root_dir, args.project_name),
            overwrite=args.overwrite,
        )

    elif args.command == "recipe":

        if args.list_only:
            import pkgutil
            import ngs_toolkit.recipes

            n = pkgutil.iter_modules(ngs_toolkit.recipes.__path__)
            print("Available ngs_toolkit recipes: '{}'.".format("', '".join([x[1] for x in n])))
        else:
            _LOGGER.info("Running recipe '{}'.".format(args.recipe_name))
            run_recipe(recipe_name=args.recipe_name, project_config=args.project_config)

    _LOGGER.debug("Completed.")


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
