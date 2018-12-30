#!/usr/bin/env python

import pytest
from .data_generator import generate_project
import os
import yaml
from peppy import Project
from ngs_toolkit.atacseq import ATACSeqAnalysis
from ngs_toolkit.general import differential_analysis, plot_differential


@pytest.fixture
def analysis(tmp_path):
    import rpy2

    # Let's make several "reallish" test projects
    to_test = list()
    project_prefix_name = "test-project"
    data_type = "ATAC-seq"
    genome_assemblies = [("human", "hg19"), ("human", "hg38"), ("mouse", "mm10")]

    n_factors = 1
    n_variables = 10000
    n_replicates = 10
    for organism, genome_assembly in [genome_assemblies[0]]:
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
        a.set_attributes()
        a.load_data()

        a.normalize(method="total")
        a.normalize(method="quantile")
        a.annotate_with_sample_metadata(quant_matrix="coverage_qnorm")

        # Sometimes the DESeq2 call fails, but it often works if repeated
        i = 1
        max_attempts = 10
        while not hasattr(a, "differential_results"):
            try:
                print("Attempt to run differential_analysis number: {}".format(i))
                differential_analysis(a)
            except rpy2.rinterface.RRuntimeError:
                i += 1
                if i == max_attempts:
                    raise RuntimeError("Tried to run differential_analysis a maximum of {} times."
                                       .format(max_attempts))

        to_test.append(a)
    return to_test[0]


@pytest.fixture
def outputs(analysis):
    prefix = os.path.join(analysis.results_dir,
                          "differential_analysis_ATAC-seq", "differential_analysis.")
    outputs = [
        prefix + "diff_region.samples.clustermap.corr.svg",
        prefix + "diff_region.samples.clustermap.svg",
        prefix + "diff_region.samples.clustermap.z0.svg",
        prefix + "diff_region.samples.sorted.clustermap.svg",
        prefix + "diff_region.samples.sorted.clustermap.z0.svg",
        prefix + "log2FoldChange.distribution.per_comparison.svg",
        prefix + "log2FoldChange.distribution.svg",
        prefix + "ma_plots.svg",
        prefix + "number_differential.directional.svg",
        prefix + "padj.distribution.per_comparison.svg",
        prefix + "padj.distribution.svg",
        prefix + "pvalue.distribution.per_comparison.svg",
        prefix + "pvalue.distribution.svg",
        prefix + "scatter_plots.svg",
        prefix + "volcano_plots.svg"]
    return outputs


class Test_differential_analysis:
    def test_no_arguments(self, analysis, outputs):
        plot_differential(analysis)
        for output in outputs:
            assert os.path.exists(output)
            assert os.stat(output).st_size > 0
