#!/usr/bin/env python

import pytest
from .data_generator import generate_project
import os
import yaml
from peppy import Project
from ngs_toolkit.atacseq import ATACSeqAnalysis
from ngs_toolkit import _CONFIG


@pytest.fixture
def analysis(tmp_path):
    # import rpy2
    tmp_path = os.path.abspath(str(tmp_path))  # for Python2

    # Let's make a "reallish" test project
    project_prefix_name = "test-project"
    data_type = "ATAC-seq"
    genome_assembly = "hg19"

    n_factors = 1
    n_variables = 10000
    n_replicates = 10

    project_name = "{}_{}_{}_{}_{}_{}".format(
        project_prefix_name, data_type, genome_assembly,
        n_factors, n_variables, n_replicates)

    generate_project(
        output_dir=tmp_path,
        project_name=project_name, genome_assembly=genome_assembly, data_type=data_type,
        n_factors=n_factors, n_replicates=n_replicates, n_variables=n_variables,
        group_fold_differences=[20])

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
    a = ATACSeqAnalysis(
        name=project_name,
        prj=Project(config),
        results_dir=os.path.join(prj_path, "results"))
    a.set_project_attributes()
    a.load_data()

    a.get_peak_gene_annotation()
    a.annotate(quant_matrix="coverage")
    a.differential_analysis(filter_support=False)

    _CONFIG['resources']['enrichr']['gene_set_libraries'] = ["GO_Biological_Process_2015"]
    a.differential_enrichment(steps=['enrichr'])

    return a


@pytest.fixture
def outputs(analysis):
    gene_set_libraries = _CONFIG['resources']['enrichr']['gene_set_libraries']
    prefix = os.path.join(analysis.results_dir,
                          "differential_analysis_ATAC-seq", "enrichments",
                          "differential_analysis.enrichr.")
    outputs = list()
    for g in gene_set_libraries:
        outputs += [
            prefix + "{}.barplot.top_5.svg".format(g),
            prefix + "{}.cluster_specific.Row_z_score.svg".format(g),
            prefix + "{}.cluster_specific.svg".format(g),
            prefix + "{}.correlation.svg".format(g),
            prefix + "{}.zscore_vs_pvalue.scatterplot.svg".format(g)]
    return outputs


# @pytest.mark.skip(reason="no way of currently testing this")
class Test_plot_differential_enrichment:
    def test_no_arguments(self, analysis, outputs):
        analysis.plot_differential_enrichment()
        for output in outputs:
            assert os.path.exists(output)
            assert os.stat(output).st_size > 0
