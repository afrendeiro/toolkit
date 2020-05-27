#!/usr/bin/env python


import os

import pandas as pd
import pytest

from .conftest import file_exists_and_not_empty
from ngs_toolkit.utils import get_this_file_or_timestamped


@pytest.fixture
def a(atac_analysis_with_input_files):
    return atac_analysis_with_input_files


def test_has_bam_files(a):
    v = [file_exists_and_not_empty(s.aligned_filtered_bam) for s in a.samples]
    assert all(v)
    v = [file_exists_and_not_empty(s.aligned_filtered_bam + ".bai") for s in a.samples]
    assert all(v)


def test_has_peak_files(a):
    v = [file_exists_and_not_empty(s.peaks) for s in a.samples]
    assert all(v)


def test_has_summit_files(a):
    v = [file_exists_and_not_empty(s.summits) for s in a.samples]
    assert all(v)


class Test_measure_coverage:
    def test_no_arguments(self, a):
        mn = get_this_file_or_timestamped(os.path.join(a.results_dir, a.name + ".matrix_raw.csv"))

        os.remove(mn)

        a.measure_coverage()

        assert file_exists_and_not_empty(mn)

    def test_distributed(self, a):
        mn = get_this_file_or_timestamped(os.path.join(a.results_dir, a.name + ".matrix_raw.csv"))

        os.remove(mn)

        a.measure_coverage(distributed=True, computing_configuration="localhost")

        # Check job files for each sample exist
        fs = list()
        for s in a.samples:
            f = os.path.join(s.sample_root, "coverage", s.name + ".peak_set_coverage.")
            for end in ["sh", "bed"]:
                fs.append(f + end)
        assert all([file_exists_and_not_empty(f) for f in fs])

        # # has to be done separately for log files because they'll empty
        # # just check for existence
        fs = list()
        for s in a.samples:
            f = os.path.join(s.sample_root, "coverage", s.name + ".peak_set_coverage.")
            for end in ["log"]:
                fs.append(f + end)
        assert all([os.path.exists(f) for f in fs])

    def test_few_samples(self, a):
        mn = get_this_file_or_timestamped(os.path.join(a.results_dir, a.name + ".matrix_raw.csv"))

        os.remove(mn)

        a.measure_coverage(samples=a.samples[:2])

        mn = get_this_file_or_timestamped(mn)
        assert file_exists_and_not_empty(mn)
        assert pd.read_csv(mn, index_col=0).shape[1] == 2

    def test_one_sample(self, a):
        mn = get_this_file_or_timestamped(os.path.join(a.results_dir, a.name + ".matrix_raw.csv"))

        os.remove(mn)

        a.measure_coverage(samples=a.samples[:1])

        mn = get_this_file_or_timestamped(mn)
        assert file_exists_and_not_empty(mn)
        assert pd.read_csv(mn, index_col=0).shape[1] == 1

    def test_missing_input_no_permissive(self, a):
        mn = get_this_file_or_timestamped(os.path.join(a.results_dir, a.name + ".matrix_raw.csv"))

        os.remove(mn)
        os.remove(a.samples[0].aligned_filtered_bam)

        with pytest.raises(IOError):
            a.measure_coverage(samples=a.samples[:1])

    def test_missing_input_all_samples(self, a):
        mn = get_this_file_or_timestamped(os.path.join(a.results_dir, a.name + ".matrix_raw.csv"))

        os.remove(mn)
        for s in a.samples:
            os.remove(s.aligned_filtered_bam)

        with pytest.raises(IOError):
            a.measure_coverage()

    def test_missing_input_with_permissive(self, a):
        mn = get_this_file_or_timestamped(os.path.join(a.results_dir, a.name + ".matrix_raw.csv"))

        os.remove(mn)
        os.remove(a.samples[0].aligned_filtered_bam)

        a.measure_coverage(samples=a.samples[:2], permissive=True)

        mn = get_this_file_or_timestamped(mn)
        assert file_exists_and_not_empty(mn)
        assert pd.read_csv(mn, index_col=0).shape[1] == 1
