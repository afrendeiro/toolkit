#!/usr/bin/env python


import os

import pandas as pd
import pytest


travis = "TRAVIS" in os.environ


@pytest.fixture
def a(atac_analysis_with_input_files):
    return atac_analysis_with_input_files


def test_has_bam_files(a):
    v = [
        os.path.exists(s.aligned_filtered_bam)
        for s in a.samples]
    assert all(v)
    v = [
        os.stat(s.aligned_filtered_bam).st_size > 0
        for s in a.samples]
    assert all(v)
    v = [
        os.path.exists(s.aligned_filtered_bam + ".bai")
        for s in a.samples]
    assert all(v)
    v = [
        os.stat(s.aligned_filtered_bam + ".bai").st_size > 0
        for s in a.samples]
    assert all(v)


def test_has_peak_files(a):
    v = [
        os.path.exists(s.peaks)
        for s in a.samples]
    assert all(v)
    v = [
        os.stat(s.peaks).st_size > 0
        for s in a.samples]
    assert all(v)


def test_has_summit_files(a):
    v = [
        os.path.exists(s.summits)
        for s in a.samples]
    assert all(v)
    v = [
        os.stat(s.summits).st_size > 0
        for s in a.samples]
    assert all(v)


class Test_measure_coverage:

    def test_no_arguments(self, a):
        mn = os.path.join(a.results_dir, a.name + ".matrix_raw.csv")

        os.remove(mn)

        a.measure_coverage()

        assert os.path.exists(mn)
        assert os.stat(mn).st_size > 0

    def test_few_samples(self, a):
        mn = os.path.join(a.results_dir, a.name + ".matrix_raw.csv")

        os.remove(mn)

        a.measure_coverage(samples=a.samples[:2])

        assert os.path.exists(mn)
        assert os.stat(mn).st_size > 0
        assert pd.read_csv(mn, index_col=0).shape[1] == 2

    def test_one_sample(self, a):
        mn = os.path.join(a.results_dir, a.name + ".matrix_raw.csv")

        os.remove(mn)

        a.measure_coverage(samples=a.samples[:1])

        assert os.path.exists(mn)
        assert os.stat(mn).st_size > 0
        assert pd.read_csv(mn, index_col=0).shape[1] == 1

    def test_missing_input_no_permissive(self, a):
        mn = os.path.join(a.results_dir, a.name + ".matrix_raw.csv")

        os.remove(mn)
        os.remove(a.samples[0].aligned_filtered_bam)

        with pytest.raises(IOError):
            a.measure_coverage(samples=a.samples[:1])

    def test_missing_input_all_samples(self, a):
        mn = os.path.join(a.results_dir, a.name + ".matrix_raw.csv")

        os.remove(mn)
        for s in a.samples:
            os.remove(s.aligned_filtered_bam)

        with pytest.raises(IOError):
            a.measure_coverage()

    def test_missing_input_with_permissive(self, a):
        mn = os.path.join(a.results_dir, a.name + ".matrix_raw.csv")

        os.remove(mn)
        os.remove(a.samples[0].aligned_filtered_bam)

        a.measure_coverage(samples=a.samples[:2], permissive=True)

        assert os.path.exists(mn)
        assert os.stat(mn).st_size > 0
        assert pd.read_csv(mn, index_col=0).shape[1] == 1
