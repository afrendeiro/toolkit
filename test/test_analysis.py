#!/usr/bin/env python


from .data_generator import generate_project


def test_analysis_creation(tmp_path):
    from ngs_toolkit.general import Analysis
    import os
    from peppy import Project
    import yaml
    import shutil

    name = "test_analysis"

    a = Analysis(name=name)
    assert a.__repr__() == "Analysis object named '{}'.".format(name)
    assert "samples" not in a.__repr__()

    # Let's make several "reallish" test projects
    project_prefix_name = "test-project"
    data_types = ["ATAC-seq", "RNA-seq", "ChIP-seq"]  # "CNV"
    genome_assemblies = [("human", "hg19"), ("human", "hg38"), ("mouse", "mm10")]

    params = {
        "ATAC-seq": {
            "n_factors": [1, 2, 3],
            "n_variables": [100, 1000, 10000],
            "n_replicates": [1, 2, 5],
        },
        "ChIP-seq": {
            "n_factors": [1, 2, 3],
            "n_variables": [100, 1000, 10000],
            "n_replicates": [1, 2, 5],
        },
        "RNA-seq": {
            "n_factors": [1, 2, 3],
            "n_variables": [100, 1000, 25000],
            "n_replicates": [1, 2, 5],
        },
    }

    for data_type in data_types:
        n_factors = params[data_type]['n_factors'][0]
        n_variables = params[data_type]['n_variables'][0]
        n_replicates = params[data_type]['n_replicates'][0]
        for organism, genome_assembly in genome_assemblies:

            project_name = "{}_{}_{}_{}_{}_{}".format(
                project_prefix_name, data_type, genome_assembly,
                n_factors, n_variables, n_replicates
            )
            n_samples = (n_factors * n_replicates) + n_factors

            generate_project(
                output_dir=tmp_path,
                project_name=project_name, genome_assembly=genome_assembly, data_type=data_type,
                n_factors=n_factors, n_replicates=n_replicates, n_variables=n_variables)

            # first edit the defaul path to the annotation sheet
            config = os.path.join(
                tmp_path, project_name, "metadata", "project_config.yaml")
            c = yaml.safe_load(open(config, 'r'))
            c['metadata']['sample_annotation'] = os.path.abspath(
                os.path.join(tmp_path, project_name, "metadata", "annotation.csv"))
            c['metadata']['comparison_table'] = os.path.abspath(
                os.path.join(tmp_path, project_name, "metadata", "comparison_table.csv"))
            yaml.safe_dump(c, open(config, "w"))

            # project and associated analysis
            prj = Project(config)
            a = Analysis(name=project_name, prj=prj)
            assert a.__repr__() == (
                "Analysis object named '{}' with {} samples of genome '{}'."
                .format(project_name, n_samples, genome_assembly))
            assert len(prj.samples) == len(a.samples)
            assert all([x == y for x, y in zip(prj.samples, a.samples)])

            shutil.rmtree(tmp_path)


def test_analysis_serialization(tmp_path):
    from ngs_toolkit.general import Analysis
    import shutil
    import os
    import numpy as np

    pickle_file = os.path.join(tmp_path, "pickle")
    a = Analysis(pickle_file=pickle_file)
    assert not os.path.exists(pickle_file)
    a.to_pickle()
    assert os.path.exists(pickle_file)
    assert os.stat(pickle_file).st_size > 0

    previous_size = os.stat(pickle_file).st_size
    a.random = np.random.random((100, 100))
    a.to_pickle()
    new_size = os.stat(pickle_file).st_size
    assert new_size > previous_size

    shutil.rmtree(tmp_path)


def test_analysis_loading(tmp_path):
    from ngs_toolkit.general import Analysis
    import shutil
    import os

    pickle_file = os.path.join(tmp_path, "pickle")
    a = Analysis(pickle_file=pickle_file)
    a.secret = "I've existed before"
    a.to_pickle()

    a2 = Analysis(pickle_file=pickle_file, from_pickle=True)
    assert a2.secret == "I've existed before"

    a3 = Analysis()
    a3.update(pickle_file)
    assert a3.secret == "I've existed before"

    a4 = Analysis(pickle_file=pickle_file).from_pickle()
    assert a4.secret == "I've existed before"

    shutil.rmtree(tmp_path)
