#!/usr/bin/env python


import pytest

from .conftest import file_exists, file_exists_and_not_empty  # , CI, RPY2


def test_call_peaks_from_comparisons(chipseq_analysis):
    chipseq_analysis.call_peaks_from_comparisons()

    for name, comp in chipseq_analysis.comparisons.items():
        files = [
            comp['prefix'] + ".homer.log",
            comp['prefix'] + ".homer.sh",
            comp['prefix'] + ".macs2.log",
            comp['prefix'] + ".macs2.sh",
        ]
        for f in files:
            assert file_exists(f)


def test_filter_peaks(chipseq_analysis_with_peaks):
    chipseq_analysis_with_peaks.filter_peaks()

    for name, comp in chipseq_analysis_with_peaks.comparisons.items():
        files = [
            comp['prefix'] + "_homer_peaks.factor.filtered.bed",
            comp['prefix'] + "_homer_peaks.factor.narrowPeak",
            comp['prefix'] + "_homer_peaks.histone.filtered.bed",
            comp['prefix'] + "_homer_peaks.histone.narrowPeak",
            comp['prefix'] + "_peaks.filtered.bed",
            comp['prefix'] + "_peaks.narrowPeak",
        ]
        for f in files:
            assert file_exists_and_not_empty(f)


@pytest.mark.xfail
def test_summarize_peaks_from_comparisons(chipseq_analysis_with_peaks):
    chipseq_analysis_with_peaks.test_summarize_peaks_from_comparisons()
    assert False


@pytest.mark.xfail
def test_get_consensus_sites(chipseq_analysis_with_peaks):
    chipseq_analysis_with_peaks.test_get_consensus_sites()
    assert False


@pytest.mark.xfail
def test_calculate_peak_support(chipseq_analysis_with_peaks):
    chipseq_analysis_with_peaks.test_calculate_peak_support()
    assert False


@pytest.mark.xfail
def test_get_supported_peaks(chipseq_analysis_with_peaks):
    chipseq_analysis_with_peaks.test_get_supported_peaks()
    assert False


@pytest.mark.xfail
def test_normalize_by_background(chipseq_analysis_with_peaks):
    chipseq_analysis_with_peaks.test_normalize_by_background()
    assert False
