#!/usr/bin/env python


import os

from .data_generator import generate_project
from ngs_toolkit.rnaseq import RNASeqAnalysis
import numpy as np
import pandas as pd
from peppy import Project
import pytest


travis = "TRAVIS" in os.environ


@pytest.fixture
def various_analysis(tmp_path):
    tmp_path = str(tmp_path)  # for Python2

    # Let's make several "reallish" test projects
    to_test = list()
    project_prefix_name = "test-project"
    data_type = "RNA-seq"
    genome_assemblies = [("human", "hg19"), ("human", "hg38"), ("mouse", "mm10")]

    n_factors = [1, 2, 3][1]
    n_variables = [100, 1000, 10000][0]
    n_replicates = [1, 2, 5][1]
    for organism, genome_assembly in genome_assemblies:
        project_name = "{}_{}_{}_{}_{}_{}".format(
            project_prefix_name,
            data_type,
            genome_assembly,
            n_factors,
            n_variables,
            n_replicates,
        )

        generate_project(
            output_dir=tmp_path,
            project_name=project_name,
            organism=organism,
            genome_assembly=genome_assembly,
            data_type=data_type,
            n_factors=n_factors,
            n_replicates=n_replicates,
            n_variables=n_variables,
        )

        # first edit the defaul path to the annotation sheet
        config = os.path.join(tmp_path, project_name, "metadata", "project_config.yaml")
        prj_path = os.path.join(tmp_path, project_name)
        os.chdir(prj_path)

        # project and associated analysis
        analysis = RNASeqAnalysis(
            name=project_name,
            prj=Project(config),
            results_dir=os.path.join(prj_path, "results"),
        )
        analysis.load_data()

        to_test.append(analysis)
    return to_test


@pytest.fixture
def analysis(tmp_path):
    tmp_path = str(tmp_path)  # for Python2

    # Let's make several "reallish" test projects
    project_prefix_name = "test-project"
    data_type = "RNA-seq"
    organism, genome_assembly = ("human", "hg38")

    n_factors = [1, 2, 3][1]
    n_variables = [100, 1000, 10000][0]
    n_replicates = [1, 2, 5][1]
    project_name = "{}_{}_{}_{}_{}_{}".format(
        project_prefix_name,
        data_type,
        genome_assembly,
        n_factors,
        n_variables,
        n_replicates,
    )

    generate_project(
        output_dir=tmp_path,
        project_name=project_name,
        organism=organism,
        genome_assembly=genome_assembly,
        data_type=data_type,
        n_factors=n_factors,
        n_replicates=n_replicates,
        n_variables=n_variables,
    )

    # first edit the defaul path to the annotation sheet
    config = os.path.join(tmp_path, project_name, "metadata", "project_config.yaml")
    # project and associated analysis
    analysis = RNASeqAnalysis(
        from_pep=config,
    )
    analysis.load_data()

    return analysis


def test_rpm_normalization(various_analysis):
    for analysis in various_analysis:
        qnorm = analysis.normalize_rpm(save=False)
        assert qnorm.dtypes.all() == np.float
        assert hasattr(analysis, "matrix_norm")
        rpm_file = os.path.join(
            analysis.results_dir, analysis.name + ".matrix_norm.csv"
        )
        assert not os.path.exists(rpm_file)
        qnorm = analysis.normalize_rpm(save=True)
        assert os.path.exists(rpm_file)
        assert os.stat(rpm_file).st_size > 0
        assert hasattr(analysis, "matrix_norm")


def test_quantile_normalization(various_analysis):
    for analysis in various_analysis:
        f = os.path.join(analysis.results_dir, analysis.name + ".matrix_norm.csv")
        qnorm_p = analysis.normalize_quantiles(implementation="Python", save=True)
        assert isinstance(qnorm_p, pd.DataFrame)
        assert hasattr(analysis, "matrix_norm")
        assert os.path.exists(f)
        assert os.stat(f).st_size > 0
        del analysis.matrix_norm
        os.remove(f)

        qnorm_r = analysis.normalize_quantiles(implementation="R", save=True)
        assert isinstance(qnorm_r, pd.DataFrame)
        assert hasattr(analysis, "matrix_norm")
        assert os.path.exists(f)
        assert os.stat(f).st_size > 0


def test_normalize(analysis):
    qnorm = analysis.normalize_rpm(save=False)
    assert isinstance(qnorm, pd.DataFrame)
    assert hasattr(analysis, "matrix_norm")
    del analysis.matrix_norm

    qnorm_d = analysis.normalize(method="rpm", save=False)
    assert isinstance(qnorm_d, pd.DataFrame)
    assert hasattr(analysis, "matrix_norm")
    assert np.array_equal(qnorm_d, qnorm)
    del analysis.matrix_norm

    qnorm = analysis.normalize_quantiles(save=False)
    assert hasattr(analysis, "matrix_norm")
    del analysis.matrix_norm

    qnorm_d = analysis.normalize(method="quantile", save=False)
    assert isinstance(qnorm_d, pd.DataFrame)
    assert hasattr(analysis, "matrix_norm")
    assert np.array_equal(qnorm_d, qnorm)
    del analysis.matrix_norm


def test_annotate_features(analysis):
    analysis.get_matrix_stats(matrix="matrix_raw")
    analysis.annotate_features(matrix="matrix_raw")
    f = os.path.join(analysis.results_dir, analysis.name + ".matrix_features.csv")
    assert hasattr(analysis, "matrix_features")
    assert os.path.exists(f)
    assert os.stat(f).st_size > 0

    cols = [
        "mean",
        "variance",
        "std_deviation",
        "dispersion",
        "qv2",
        "amplitude",
        "iqr",
    ]  # from stats

    assert all([c in analysis.matrix_features.columns.tolist() for c in cols])


def test_plot_expression_characteristics(various_analysis):
    for analysis in various_analysis:
        analysis.normalize()
        analysis.plot_expression_characteristics()
        assert os.path.exists(os.path.join(analysis.results_dir, "quality_control"))
