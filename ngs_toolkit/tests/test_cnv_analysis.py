#!/usr/bin/env python

import os

import pytest
import pandas as pd

from .conftest import file_exists_and_not_empty, STAP, DNACOPY  # , CI, RPY2


@pytest.mark.xfail
def test__copy_cnv_profile_plots(cnv_analysis):
    with cnv_analysis as an:
        an._copy_cnv_profile_plots()
    assert False


def test_get_cnv_data(cnv_analysis_with_inputs):
    with cnv_analysis_with_inputs as an:
        an.get_cnv_data()

        p = os.path.join(an.results_dir, an.name)
        files = [
            p + ".10kb.matrix_raw.csv",
            p + ".100kb.matrix_raw.csv",
            p + ".1000kb.matrix_raw.csv"]
        for f in files:
            assert file_exists_and_not_empty(f)
            assert pd.read_csv(f, index_col=0).sum().sum() == 0


def test_normalize(cnv_analysis):
    with cnv_analysis as an:
        an.normalize()

        p = os.path.join(an.results_dir, an.name)
        files = [
            p + ".10kb.matrix_norm.csv",
            p + ".100kb.matrix_norm.csv",
            p + ".1000kb.matrix_norm.csv"]
        for f in files:
            assert file_exists_and_not_empty(f)
            assert pd.read_csv(f, index_col=0).sum().sum() != 0


def test_plot_all_data(cnv_analysis):
    with cnv_analysis as an:
        an.plot_all_data(matrix='matrix_raw')

        for res in an.resolutions:
            p = os.path.join(an.results_dir, an.name + "." + res + ".all_data.full_data")
            files = [
                p + ".fillna.clustermap.svg",
                p + ".heatmap.svg"]
            for f in files:
                assert file_exists_and_not_empty(f)


def test_plot_stats_per_chromosome(cnv_analysis):
    with cnv_analysis as an:
        an.plot_stats_per_chromosome(matrix="matrix_raw")

        for res in an.resolutions:
            for t in ['mean', 'variation']:
                p = os.path.join(an.results_dir, an.name + "." + res + ".all_data." + t + "_per_chrom")
                files = [
                    p + ".no_sex_chroms.zscore.svg",
                    p + ".no_sex_chroms.svg",
                    p + ".svg"]
                for f in files:
                    assert file_exists_and_not_empty(f)


# @pytest.mark.skipif(not STAP or not DNACOPY, reason="STAP and DNACopy R libraries are required to perform segmentation.")
@pytest.mark.xfail
def test_segment_genome(cnv_analysis):
    cnv_analysis.segment_genome()
    assert False


@pytest.mark.xfail
def test_annotate_with_chrom_bands(cnv_analysis):
    cnv_analysis.annotate_with_chrom_bands()
    assert False


@pytest.mark.xfail
def test_plot_segmentation_stats(cnv_analysis):
    cnv_analysis.plot_segmentation_stats()
    assert False
