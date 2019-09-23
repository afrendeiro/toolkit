#!/usr/bin/env python


class Test_annotate_samples:
    def test_no_arguments(self, analysis_normalized):
        analysis_normalized.annotate_samples()

    def test_matrix_raw(self, atac_analysis):
        atac_analysis.annotate_samples(matrix="matrix_raw")


def test_get_matrix(atac_analysis):
    import numpy as np
    import pandas as pd

    matrix = atac_analysis.get_matrix(matrix="matrix_raw")
    assert np.array_equal(matrix.values, atac_analysis.matrix_raw.values)
    assert (matrix == atac_analysis.matrix_raw).all().all()
    atac_analysis.dummy = atac_analysis.matrix_raw + 1
    matrix = atac_analysis.get_matrix(matrix="dummy")
    assert (matrix == (atac_analysis.matrix_raw + 1)).all().all()

    matrix = atac_analysis.get_matrix(matrix=atac_analysis.matrix_raw)
    assert np.array_equal(matrix.values, atac_analysis.matrix_raw.values)
    assert (matrix == atac_analysis.matrix_raw).all().all()
    atac_analysis.dummy = atac_analysis.matrix_raw + 1
    matrix = atac_analysis.get_matrix(matrix="dummy")
    assert (matrix == (atac_analysis.matrix_raw + 1)).all().all()

    # sample subssetting
    matrix = atac_analysis.get_matrix(
        matrix="matrix_raw", samples=atac_analysis.samples[:2])
    assert (pd.Series([
        s.name
        for s in atac_analysis.samples[:2]]) == matrix.columns).all()
