#!/usr/bin/env python

import pytest
import os
import yaml
from .test_analysis import generate_project
from peppy import Project
from ngs_toolkit.atacseq import ATACSeqAnalysis
import shutil
import pybedtools
import pandas as pd
import numpy as np


@pytest.fixture
def get_test_analysis(tmp_path):
    # Let's make several "reallish" test projects
    to_test = list()
    project_prefix_name = "test-project"
    data_type = "ATAC-seq"
    genome_assemblies = [("human", "hg19"), ("human", "hg38"), ("mouse", "mm10")]

    n_factors = [1, 2, 3][0]
    n_variables = [100, 1000, 10000][0]
    n_replicates = [1, 2, 5][0]
    for organism, genome_assembly in genome_assemblies:
        project_name = "{}_{}_{}_{}_{}_{}".format(
            project_prefix_name, data_type, genome_assembly,
            n_factors, n_variables, n_replicates)

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

        prj_path = os.path.join(tmp_path, project_name)
        os.chdir(prj_path)

        # project and associated analysis
        analysis = ATACSeqAnalysis(
            name=project_name,
            prj=Project(config),
            results_dir=os.path.join(prj_path, "results"))
        analysis.load_data()

        to_test.append(analysis)
    return to_test


def test_consensus_set_loading(get_test_analysis):
    for analysis in get_test_analysis:
        # load up consesnsus set
        assert hasattr(analysis, "sites")
        assert isinstance(analysis.sites, pybedtools.BedTool)


def test_coverage_matrix_loading(get_test_analysis):
    for analysis in get_test_analysis:
        # load up coverage matrix
        assert hasattr(analysis, "coverage")
        assert isinstance(analysis.coverage, pd.DataFrame)
        assert analysis.coverage.dtypes.all() == int


def test_setting_consensus_set(get_test_analysis):
    for analysis in get_test_analysis:
        # setting a new consensus set
        peaks = os.path.join(analysis.results_dir, analysis.name + "_peak_set.bed")
        analysis.set_consensus_sites(peaks)
        assert hasattr(analysis, "sites")
        sites = pd.read_csv(peaks, header=None)
        assert len(analysis.sites) == sites.shape[0]


def test_get_matrix(get_test_analysis):
    for analysis in get_test_analysis:
        matrix = analysis.get_matrix()
        assert np.array_equal(matrix.values, analysis.coverage.values)
        assert (matrix == analysis.coverage).all().all()
        analysis.dummy = analysis.coverage + 1
        matrix = analysis.get_matrix(matrix_name="dummy")
        assert (matrix == (analysis.coverage + 1)).all().all()


def test_rpm_normalization(get_test_analysis):
    for analysis in get_test_analysis:
        qnorm = analysis.normalize_coverage_rpm(save=False)
        assert qnorm.dtypes.all() == np.float
        rpm_file = os.path.join(analysis.results_dir, analysis.name + "_peaks.coverage_rpm.csv")
        assert not os.path.exists(rpm_file)
        qnorm = analysis.normalize_coverage_rpm(save=True)
        assert os.path.exists(rpm_file)
        assert os.stat(rpm_file).st_size > 0


def test_quantile_normalization(get_test_analysis):
    for analysis in get_test_analysis:
        qnorm_p = analysis.normalize_coverage_quantiles(implementation="Python", save=False)
        qnorm_r = analysis.normalize_coverage_quantiles(implementation="R", save=False)

        import scipy
        cors = list()
        for col in qnorm_p.columns:
            cors.append(scipy.stats.pearsonr(qnorm_p[col], qnorm_r[col])[0])
        assert all(np.array(cors) > 0.99)


# TODO: test cqn normalization


def test_normalize(get_test_analysis):
    for analysis in get_test_analysis:
        qnorm = analysis.normalize_coverage_rpm(save=False)
        qnorm_d = analysis.normalize(method="total", save=False)
        assert np.array_equal(qnorm_d, qnorm)
        qnorm = analysis.normalize_coverage_quantiles(save=False)
        qnorm_d = analysis.normalize(method="quantile", save=False)
        assert np.array_equal(qnorm_d, qnorm)
        # TODO: add cqn normalization


# TODO: test region set annotation


def test_plot_raw_coverage(get_test_analysis):
    for analysis in get_test_analysis:
        analysis.plot_raw_coverage()
        output = os.path.join(analysis.results_dir, analysis.name + ".raw_counts.violinplot.svg")
        assert os.path.exists(output)
        assert os.stat(output).st_size > 0

        attr = "a"
        analysis.plot_raw_coverage(by_attribute=attr)
        output = os.path.join(analysis.results_dir, analysis.name + ".raw_counts.violinplot.by_{}.svg".format(attr))
        assert os.path.exists(output)
        assert os.stat(output).st_size > 0


# def test_cleanup(get_test_analysis):
#     for analysis in get_test_analysis:
#         shutil.rmtree(analysis.results_dir)
