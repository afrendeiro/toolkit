#!/usr/bin/env python

import os

import pytest

from .conftest import file_exists, file_exists_and_not_empty


@pytest.fixture
def outputs(atac_analysis):
    output_dir = os.path.join(atac_analysis.results_dir, "differential_analysis_ATAC-seq")
    prefix = os.path.join(output_dir, "differential_analysis.")
    outputs = [
        os.path.join(output_dir, "Factor_a_2vs1"),
        os.path.join(
            output_dir,
            "Factor_a_2vs1",
            "differential_analysis.deseq_result.Factor_a_2vs1.csv",
        ),
        prefix + "comparison_table.tsv",
        prefix + "count_matrix.tsv",
        prefix + "deseq_result.all_comparisons.csv",
        prefix + "experiment_matrix.tsv",
    ]
    return outputs


# @pytest.fixture
# def outputs_no_subdirectories(analysis):
#     output_dir = os.path.join(analysis.results_dir, "differential_analysis_ATAC-seq")
#     prefix = os.path.join(output_dir, "differential_analysis.")
#     outputs = [
#         prefix + "deseq_result.Factor_a_2vs1.csv",
#         prefix + "comparison_table.tsv",
#         prefix + "count_matrix.tsv",
#         prefix + "deseq_result.all_comparisons.csv",
#         prefix + "experiment_matrix.tsv"]
#     return outputs


class Test_differential_analysis:
    def test_simple_design(self, atac_analysis, outputs):
        import pandas as pd

        atac_analysis.differential_analysis(filter_support=False)
        assert file_exists(
            os.path.join(atac_analysis.results_dir, "differential_analysis_ATAC-seq")
        )
        assert file_exists(outputs[0])
        assert os.path.isdir(outputs[0])
        for output in outputs[1:]:
            assert file_exists_and_not_empty(output)
        assert hasattr(atac_analysis, "differential_results")
        assert isinstance(atac_analysis.differential_results, pd.DataFrame)
        assert atac_analysis.differential_results.index.str.startswith("chr").all()
        assert atac_analysis.differential_results.index.name == "index"
        cols = [
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "stat",
            "pvalue",
            "padj",
            "comparison_name",
        ]
        assert atac_analysis.differential_results.columns.tolist() == cols

    def test_complex_design(self, atac_analysis, outputs):
        import pandas as pd

        atac_analysis.differential_analysis(filter_support=False)
        assert file_exists(
            os.path.join(atac_analysis.results_dir, "differential_analysis_ATAC-seq")
        )
        assert file_exists(outputs[0])
        assert os.path.isdir(outputs[0])
        for output in outputs[1:]:
            assert file_exists_and_not_empty(output)
        assert hasattr(atac_analysis, "differential_results")
        assert isinstance(atac_analysis.differential_results, pd.DataFrame)
        assert atac_analysis.differential_results.index.str.startswith("chr").all()
        assert atac_analysis.differential_results.index.name == "index"
        cols = [
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "stat",
            "pvalue",
            "padj",
            "comparison_name",
        ]
        assert atac_analysis.differential_results.columns.tolist() == cols

    # def test_no_subdirectories(self, atac_analysis, outputs):
    #     atac_analysis.differential_analysis()
    #     assert file_exists(
    #         os.path.join(atac_analysis.results_dir, "differential_analysis_ATAC-seq"))
    #     assert file_exists(outputs[0])
    #     assert os.path.isdir(outputs[0])
    #     for output in outputs[1:]:
    #         assert file_exists(output)
    #         assert os.stat(output).st_size > 0
