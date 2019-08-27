#!/usr/bin/env python

from .conftest import file_exists, file_exists_and_not_empty


def test_generate_report(atac_analysis):
    report = atac_analysis._format_string_with_attributes(
        "{root_dir}/{name}.analysis_report.html")
    assert not file_exists(report)
    atac_analysis.generate_report()
    assert file_exists_and_not_empty(report)
