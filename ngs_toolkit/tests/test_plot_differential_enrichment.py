#!/usr/bin/env python

import os

import pytest
from .conftest import file_exists_and_not_empty, R, R_REASON


@pytest.fixture
def outputs(analysis_with_differential_enrichment):
    # gene_set_libraries = _CONFIG['resources']['enrichr']['gene_set_libraries']
    gene_set_libraries = ["GO_Biological_Process_2015", "NCI-Nature_2016"]
    prefix = os.path.join(
        analysis_with_differential_enrichment.results_dir,
        "differential_analysis_ATAC-seq",
        "enrichments",
        "differential_analysis.enrichr.",
    )
    outputs = list()
    for g in gene_set_libraries:
        outputs += [
            prefix + "{}.barplot.top_5.svg".format(g),
            prefix + "{}.cluster_specific.Row_z_score.svg".format(g),
            prefix + "{}.cluster_specific.svg".format(g),
            prefix + "{}.correlation.svg".format(g),
            prefix + "{}.zscore_vs_pvalue.scatterplot.svg".format(g),
        ]
    return outputs


@pytest.mark.skipif(
    not R,
    reason=R_REASON)
class Test_plot_differential_enrichment:
    def test_no_arguments(self, analysis_with_differential_enrichment, outputs):
        analysis_with_differential_enrichment.plot_differential_enrichment()
        for output in outputs:
            assert file_exists_and_not_empty(output)
