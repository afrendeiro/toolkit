#!/usr/bin/env python

import os
import shutil

import pytest


@pytest.fixture
def unsup_outputs(atac_analysis_many_factors):
    prefix = os.path.join(
        atac_analysis_many_factors.results_dir,
        "unsupervised_analysis_{}".format(atac_analysis_many_factors.data_type),
        atac_analysis_many_factors.name + ".all_{}s.".format(atac_analysis_many_factors.var_unit_name),
    )
    outputs = [
        prefix + "isomap.svg",
        prefix + "locallylinearembedding.svg",
        prefix + "mds.svg",
        prefix + "pca.explained_variance.csv",
        prefix + "pca.explained_variance.svg",
        prefix + "pca.svg",
        prefix + "pca.variable_principle_components_association.csv",
        prefix + "pca.variable_principle_components_association.p_value.masked.svg",
        prefix + "pca.variable_principle_components_association.p_value.svg",
        prefix + "pca.variable_principle_components_association.adj_pvalue.masked.svg",
        prefix + "pca.variable_principle_components_association.adj_pvalue.svg",
        prefix + "pearson_correlation.clustermap.svg",
        prefix + "spearman_correlation.clustermap.svg",
        prefix + "spectralembedding.svg",
        prefix + "tsne.svg",
    ]
    return outputs


class TestUnsupervisedAnalysis:
    def test_no_arguments(self, atac_analysis_many_factors, unsup_outputs):
        # no arguments
        atac_analysis_many_factors.unsupervised_analysis()
        for output in unsup_outputs:
            assert os.path.exists(output)
            assert os.stat(output).st_size > 0

    def test_matrix_with_no_group_attributes(self, atac_analysis_many_factors):
        atac_analysis_many_factors.group_attributes = []
        with pytest.raises(ValueError):
            atac_analysis_many_factors.unsupervised_analysis()

    def test_matrix_with_no_multiindex(self, atac_analysis_many_factors, unsup_outputs):
        atac_analysis_many_factors.unsupervised_analysis(matrix="matrix_raw")
        for output in unsup_outputs:
            assert os.path.exists(output)
            assert os.stat(output).st_size > 0

    def test_matrix_with_no_multiindex_no_sample_attributes(self, atac_analysis):
        atac_analysis.sample_attributes = []
        with pytest.raises(ValueError):
            atac_analysis.unsupervised_analysis(matrix="matrix_raw")

    def test_matrix_with_no_multiindex2(self, atac_analysis_many_factors):
        atac_analysis_many_factors.unsupervised_analysis(matrix="matrix_raw")
        assert os.path.exists(
            os.path.join(atac_analysis_many_factors.results_dir, "unsupervised_analysis_ATAC-seq")
        )

    def test_various_matrices(self, atac_analysis_many_factors, unsup_outputs):
        for matrix in ["matrix_raw", "matrix_norm"]:
            atac_analysis_many_factors.annotate_samples(matrix=matrix)
            atac_analysis_many_factors.unsupervised_analysis()
            for output in unsup_outputs:
                assert os.path.exists(output)
                assert os.stat(output).st_size > 0
            shutil.rmtree(
                os.path.join(atac_analysis_many_factors.results_dir, "unsupervised_analysis_ATAC-seq")
            )
        # analysis.annotate_samples(matrix="coverage_qnorm")

    def test_too_low_numbers_of_samples_error(self, atac_analysis_many_factors):
        for i in range(2):
            with pytest.raises(ValueError):
                atac_analysis_many_factors.unsupervised_analysis(samples=atac_analysis_many_factors.samples[:i])
            assert os.path.exists(
                os.path.join(atac_analysis_many_factors.results_dir, "unsupervised_analysis_ATAC-seq")
            )
            shutil.rmtree(
                os.path.join(atac_analysis_many_factors.results_dir, "unsupervised_analysis_ATAC-seq")
            )

    def test_low_samples_no_manifolds(self, atac_analysis_many_factors):
        prefix = os.path.join(
            atac_analysis_many_factors.results_dir,
            "unsupervised_analysis_ATAC-seq",
            atac_analysis_many_factors.name + ".all_{}s.".format(atac_analysis_many_factors.var_unit_name),
        )
        outputs2 = [
            prefix + "mds.svg",
            prefix + "pca.explained_variance.csv",
            prefix + "pca.explained_variance.svg",
            prefix + "pca.svg",
            prefix + "pearson_correlation.clustermap.svg",
            prefix + "spearman_correlation.clustermap.svg",
            prefix + "tsne.svg",
            prefix + "pca.variable_principle_components_association.csv",
        ]
        not_outputs = [
            prefix + "isomap.svg",
            prefix + "locallylinearembedding.svg",
            prefix + "spectralembedding.svg",
            prefix + "pca.variable_principle_components_association.p_value.masked.svg",
            prefix
            + "pca.variable_principle_components_association.adj_pvalue.masked.svg",
            prefix + "pca.variable_principle_components_association.p_value.svg",
            prefix + "pca.variable_principle_components_association.adj_pvalue.svg",
        ]
        # here I'm picking the first and last samples just to make sure
        # they are from different values of attributes `a` and `b`
        atac_analysis_many_factors.unsupervised_analysis(
            samples=[atac_analysis_many_factors.samples[0]] + [atac_analysis_many_factors.samples[-1]])
        for output in outputs2:
            assert os.path.exists(output)
            assert os.stat(output).st_size > 0
        for output in not_outputs:
            assert not os.path.exists(output)

    # def test_high_samples_varying_all_outputs(self, atac_analysis_many_factors, outputs):
    #     for i in range(4, len(atac_analysis_many_factors.samples), 2):
    #         print(i)
    #         atac_analysis_many_factors.unsupervised_analysis(samples=atac_analysis_many_factors.samples[i:])
    #         for output in outputs:
    #             assert os.path.exists(output)
    #             assert os.stat(output).st_size > 0
    #         shutil.rmtree(os.path.join(atac_analysis_many_factors.results_dir, "unsupervised_analysis_ATAC-seq"))

    def test_no_plotting_attributes(self, atac_analysis_many_factors):
        with pytest.raises(ValueError):
            atac_analysis_many_factors.unsupervised_analysis(attributes_to_plot=[])
        assert os.path.exists(
            os.path.join(atac_analysis_many_factors.results_dir, "unsupervised_analysis_ATAC-seq")
        )

    def test_various_plotting_attributes(self, analysis_annotated, unsup_outputs):
        prefix = os.path.join(
            analysis_annotated.results_dir,
            "unsupervised_analysis_ATAC-seq",
            analysis_annotated.name + ".all_{}s.".format(analysis_annotated.var_unit_name),
        )
        not_outputs = [
            prefix + "pca.variable_principle_components_association.p_value.masked.svg",
            prefix + "pca.variable_principle_components_association.p_value.svg",
            prefix
            + "pca.variable_principle_components_association.adj_pvalue.masked.svg",
            prefix + "pca.variable_principle_components_association.adj_pvalue.svg",
        ]
        for i in range(1, len(analysis_annotated.group_attributes)):
            analysis_annotated.unsupervised_analysis(
                attributes_to_plot=analysis_annotated.group_attributes[:i]
            )
            for output in unsup_outputs:
                if output not in not_outputs:
                    assert os.path.exists(output)
                    assert os.stat(output).st_size > 0
            for output in not_outputs:
                assert not os.path.exists(output)
            shutil.rmtree(
                os.path.join(analysis_annotated.results_dir, "unsupervised_analysis_ATAC-seq")
            )

    def test_various_plot_prefixes_attributes(self, atac_analysis_many_factors, unsup_outputs):
        atac_analysis_many_factors.unsupervised_analysis(plot_prefix="test")
        for output in unsup_outputs:
            assert os.path.exists(
                output.replace("all_{}s".format(atac_analysis_many_factors.var_unit_name), "test")
            )
            assert (
                os.stat(
                    output.replace("all_{}s".format(atac_analysis_many_factors.var_unit_name), "test")
                ).st_size
                > 0
            )

    def test_standardized_matrix(self, atac_analysis_many_factors, unsup_outputs):
        atac_analysis_many_factors.unsupervised_analysis(standardize_matrix=True)
        for output in unsup_outputs:
            assert os.path.exists(output)
            assert os.stat(output).st_size > 0
