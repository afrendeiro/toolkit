#!/usr/bin/env python

import pytest
from .data_generator import generate_project
import os
import yaml
from peppy import Project
from ngs_toolkit.atacseq import ATACSeqAnalysis


@pytest.fixture
def get_test_analysis(tmp_path):
    # Let's make several "reallish" test projects
    to_test = list()
    project_prefix_name = "test-project"
    data_type = "ATAC-seq"
    genome_assemblies = [("human", "hg19"), ("human", "hg38"), ("mouse", "mm10")]

    n_factors = 1
    n_variables = 100
    n_replicates = 10
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
        analysis.set_attributes()
        analysis.load_data()

        analysis.normalize(method="quantile")

        to_test.append(analysis)
    return to_test


def test_annotate_with_sample_metadata(get_test_analysis):
    for analysis in get_test_analysis:
        analysis.annotate_with_sample_metadata(quant_matrix="coverage_qnorm")


def test_unsupervised_analysis(get_test_analysis):
    from ngs_toolkit.general import unsupervised_analysis
    for analysis in get_test_analysis:
        analysis.annotate_with_sample_metadata(quant_matrix="coverage_qnorm")

        unsupervised_analysis(analysis)

        prefix = os.path.join(
            analysis.results_dir, "unsupervised_analysis_ATAC-seq", analysis.name)
        outputs = [
            prefix + ".all_sites.isomap.svg",
            prefix + ".all_sites.locallylinearembedding.svg",
            prefix + ".all_sites.mds.svg",
            prefix + ".all_sites.pca.explained_variance.csv",
            prefix + ".all_sites.pca.explained_variance.svg",
            prefix + ".all_sites.pca.svg",
            prefix + ".all_sites.pca.variable_principle_components_association.csv",
            prefix + ".all_sites.pearson_correlation.clustermap.svg",
            prefix + ".all_sites.spearman_correlation.clustermap.svg",
            prefix + ".all_sites.spectralembedding.svg",
            prefix + ".all_sites.tsne.svg"]
        for output in outputs:
            assert os.path.exists(output)
            assert os.stat(output).st_size > 0
