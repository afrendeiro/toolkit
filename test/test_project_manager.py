#!/usr/bin/env python


def test_cli_parsing():
    import pytest
    from ngs_toolkit.project_manager import parse_arguments

    with pytest.raises(SystemExit):
        parse_arguments()
        parse_arguments("")
        parse_arguments("--help")
        parse_arguments("create")
        parse_arguments("recipe")
    args = parse_arguments("create asd")
    assert args.command == "create"


def test_project_creation(tmp_path):
    from ngs_toolkit import _CONFIG
    from ngs_toolkit.project_manager import create_project
    import os
    import pandas as pd
    import shutil

    tmp_path = str(tmp_path)  # for Python2

    project_name = "test_project"
    annotation_vars = [
        "sample_name", "toggle", "pass_qc", "protocol", "library",
        "cell_line", "cell_type", "condition",
        "experimental_batch", "experiment_name", "replicate",
        "organism", "flowcell", "lane", "BSF_name", "data_source"]

    genome_assemblies = {k: v for x in _CONFIG["preferences"]["default_genome_assemblies"]
                         for k, v in x.items()}
    create_project(
        project_name, tmp_path, genome_assemblies=genome_assemblies, overwrite=True)

    expected_files = [
        os.path.join(tmp_path, project_name, ".git"),
        os.path.join(tmp_path, project_name, "metadata"),
        os.path.join(tmp_path, project_name, "metadata", "project_config.yaml"),
        os.path.join(tmp_path, project_name, "metadata", "annotation.csv"),
        os.path.join(tmp_path, project_name, "metadata", "merge_table.csv"),
        os.path.join(tmp_path, project_name, "metadata", "comparison_table.csv"),
    ]
    for f in expected_files:
        assert os.path.exists(f)

    df = pd.read_csv(os.path.join(tmp_path, project_name, "metadata", "annotation.csv"))
    assert df.shape == (0, len(annotation_vars))
    assert all(c in df.columns for c in annotation_vars)

    shutil.rmtree(tmp_path)
