#!/usr/bin/env python

import pytest
from .data_generator import generate_project
import os
import yaml
from peppy import Project
from ngs_toolkit.atacseq import ATACSeqAnalysis


@pytest.fixture
def analysis(tmp_path):
    tmp_path = str(tmp_path)  # for Python2

    # Let's make several "reallish" test projects
    to_test = list()
    project_prefix_name = "test-project"
    data_type = "ATAC-seq"
    genome_assemblies = [("human", "hg19"), ("human", "hg38"), ("mouse", "mm10")]

    n_factors = 2
    n_variables = 10000
    n_replicates = 6
    for organism, genome_assembly in [genome_assemblies[0]]:
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
        c['metadata']['output_dir'] = os.path.abspath(tmp_path)
        c['metadata']['sample_annotation'] = os.path.abspath(
            os.path.join(tmp_path, project_name, "metadata", "annotation.csv"))
        c['metadata']['comparison_table'] = os.path.abspath(
            os.path.join(tmp_path, project_name, "metadata", "comparison_table.csv"))
        yaml.safe_dump(c, open(config, "w"))

        prj_path = os.path.join(tmp_path, project_name)
        os.chdir(prj_path)

        # project and associated analysis
        a = ATACSeqAnalysis(
            name=project_name,
            prj=Project(config),
            results_dir=os.path.join(prj_path, "results"))
        a.set_project_attributes()
        a.load_data()

        a.normalize(method="rpm")
        a.annotate_with_sample_attributes()

        to_test.append(a)
    return to_test[0]


class Test_annotate_with_sample_attributes:
    def test_no_arguments(self, analysis):
        analysis.annotate_with_sample_attributes()


def test_get_matrix(analysis):
    import numpy as np
    import pandas as pd

    matrix = analysis.get_matrix(matrix="matrix_raw")
    assert np.array_equal(matrix.values, analysis.matrix_raw.values)
    assert (matrix == analysis.matrix_raw).all().all()
    analysis.dummy = analysis.matrix_raw + 1
    matrix = analysis.get_matrix(matrix="dummy")
    assert (matrix == (analysis.matrix_raw + 1)).all().all()

    matrix = analysis.get_matrix(matrix=analysis.matrix_raw)
    assert np.array_equal(matrix.values, analysis.matrix_raw.values)
    assert (matrix == analysis.matrix_raw).all().all()
    analysis.dummy = analysis.matrix_raw + 1
    matrix = analysis.get_matrix(matrix="dummy")
    assert (matrix == (analysis.matrix_raw + 1)).all().all()

    # sample subssetting
    matrix = analysis.get_matrix(matrix="matrix_raw", samples=analysis.samples[:2])
    assert (pd.Series([s.name for s in analysis.samples[:2]]) == matrix.columns).all()
