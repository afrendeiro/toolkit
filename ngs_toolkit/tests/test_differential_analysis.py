#!/usr/bin/env python

import os

import pytest

from .conftest import file_exists, file_exists_and_not_empty, R, R_REASON


# Note:
# The DESeq2 1.24.0 version in Debian archives
# differs from the DESeq2 1.24.0 version in bioconductor version 3.9
# If estimateDispersions with default fitType="parametric" fails,
# (as often happens with the quickly generated synthetic data from tests),
# it tries to use local fit using the locfit package, but in Debian
# version this is not a valid choice of fit, causing failure.
# Due to this, and since I'm using Debian packages for faster testing
# I'm manually setting fitType="mean" for testing only.


@pytest.fixture
def outputs(atac_analysis):
    output_dir = os.path.join(atac_analysis.results_dir, "differential_analysis_ATAC-seq")
    prefix = os.path.join(output_dir, "differential_analysis.")
    outs = [
        os.path.join(output_dir, "Factor_A_2vs1"),
        os.path.join(
            output_dir,
            "Factor_A_2vs1",
            "differential_analysis.deseq_result.Factor_A_2vs1.csv",
        ),
        prefix + "comparison_table.tsv",
        prefix + "count_matrix.tsv",
        prefix + "deseq_result.all_comparisons.csv",
        prefix + "experiment_matrix.tsv",
    ]
    return outs


# @pytest.fixture
# def outputs_no_subdirectories(analysis):
#     output_dir = os.path.join(analysis.results_dir, "differential_analysis_ATAC-seq")
#     prefix = os.path.join(output_dir, "differential_analysis.")
#     outputs = [
#         prefix + "deseq_result.Factor_A_2vs1.csv",
#         prefix + "comparison_table.tsv",
#         prefix + "count_matrix.tsv",
#         prefix + "deseq_result.all_comparisons.csv",
#         prefix + "experiment_matrix.tsv"]
#     return outputs


@pytest.mark.skipif(
    not R,
    reason=R_REASON)
def test_deseq_functionality():
    import pandas as pd
    from ngs_toolkit.utils import recarray2pandas_df

    from rpy2.robjects import numpy2ri, pandas2ri, r
    from rpy2.robjects.packages import importr
    numpy2ri.activate()
    pandas2ri.activate()

    importr("DESeq2")

    dds = r.makeExampleDESeqDataSet()
    dds = r.estimateSizeFactors(dds)
    dds = r.estimateDispersions(dds)
    dds = r.nbinomWaldTest(dds)
    res = recarray2pandas_df(r("as.data.frame")(r("DESeq2::results")(dds)))
    assert isinstance(res, pd.DataFrame)

    dds = r.makeExampleDESeqDataSet()
    dds = r.DESeq(dds)
    res = recarray2pandas_df(r("as.data.frame")(r("DESeq2::results")(dds)))
    assert isinstance(res, pd.DataFrame)


@pytest.mark.skipif(
    not R,
    reason=R_REASON)
class Test_differential_analysis:
    def test_simple_design(self, atac_analysis, outputs):
        import pandas as pd

        atac_analysis.differential_analysis()
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

        atac_analysis.differential_analysis()
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
