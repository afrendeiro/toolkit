#!/usr/bin/env python

import os


def test_generate_report(atac_analysis):
    report = atac_analysis._format_string_with_attributes(
        "{root_dir}/{name}.analysis_report.html")
    assert not os.path.exists(report)
    atac_analysis.generate_report()
    assert os.path.exists(report)
    assert os.stat(report).st_size > 0
