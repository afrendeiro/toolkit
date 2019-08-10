#!/usr/bin/env python

import os

import pytest


@pytest.fixture
def outputs(analysis_with_differential):
    prefix = os.path.join(
        analysis_with_differential.results_dir, "differential_analysis_ATAC-seq", "enrichments"
    )
    outputs = [
        os.path.join(
            prefix, "Factor_a_2vs1.down/differential_analysis.gene_symbols.txt"
        ),
        os.path.join(prefix, "Factor_a_2vs1.down/differential_analysis_regions.bed"),
        os.path.join(prefix, "Factor_a_2vs1.down/differential_analysis_regions.tsv"),
        os.path.join(
            prefix, "Factor_a_2vs1.down/differential_analysis_regions.enrichr.csv"
        ),
        os.path.join(prefix, "Factor_a_2vs1.up/differential_analysis.gene_symbols.txt"),
        os.path.join(prefix, "Factor_a_2vs1.up/differential_analysis_regions.bed"),
        os.path.join(prefix, "Factor_a_2vs1.up/differential_analysis_regions.tsv"),
        os.path.join(
            prefix, "Factor_a_2vs1.up/differential_analysis_regions.enrichr.csv"
        ),
        os.path.join(prefix, "differential_analysis.enrichr.csv"),
    ]
    return outputs


# @pytest.mark.skip(reason="no way of currently testing this")
class Test_differential_enrichment:
    def test_no_arguments(self, analysis_with_differential, outputs):
        analysis_with_differential.differential_enrichment(steps=["enrichr"])
        for output in outputs:
            assert os.path.exists(output)
            assert os.stat(output).st_size > 0
