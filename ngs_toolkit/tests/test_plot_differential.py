#!/usr/bin/env python

import os

import pytest
from .conftest import file_exists_and_not_empty, R, R_REASON


@pytest.fixture
def outputs(analysis_with_differential):
    prefix = os.path.join(
        analysis_with_differential.results_dir,
        "differential_analysis_ATAC-seq",
        "differential_analysis.",
    )
    outputs = [
        # prefix + "diff_region.samples.clustermap.corr.svg",
        # prefix + "diff_region.samples.clustermap.svg",
        # prefix + "diff_region.samples.clustermap.z0.svg",
        # prefix + "diff_region.samples.sorted.clustermap.svg",
        # prefix + "diff_region.samples.sorted.clustermap.z0.svg",
        prefix + "log2FoldChange.distribution.per_comparison.svg",
        prefix + "log2FoldChange.distribution.svg",
        prefix + "ma_plots.svg",
        prefix + "number_differential.directional.svg",
        prefix + "padj.distribution.per_comparison.svg",
        prefix + "padj.distribution.svg",
        prefix + "pvalue.distribution.per_comparison.svg",
        prefix + "pvalue.distribution.svg",
        prefix + "scatter_plots.svg",
        prefix + "volcano_plots.svg",
    ]
    return outputs


@pytest.mark.skipif(
    not R,
    reason=R_REASON)
class Test_plot_differential:
    def test_no_arguments(self, analysis_with_differential, outputs):
        analysis_with_differential.plot_differential()
        for output in outputs:
            assert file_exists_and_not_empty(output)
