#!/usr/bin/env python


import os

import numpy as np
import pandas as pd
import pybedtools
import pytest

from ngs_toolkit import _CONFIG
from ngs_toolkit.atacseq import ATACSeqAnalysis
from .conftest import file_exists, file_exists_and_not_empty, CI, R, R_REASON


def test_get_consensus_sites(various_analysis):
    for analysis in various_analysis:
        with pytest.raises(IOError):
            analysis.get_consensus_sites()


def test_get_supported_peaks(various_analysis):
    for analysis in various_analysis:
        analysis.support = pd.DataFrame(
            np.random.binomial(1, 0.4, size=analysis.matrix_raw.shape),
            index=analysis.matrix_raw.index,
            columns=analysis.matrix_raw.columns,
        )
        fs = analysis.get_supported_peaks(samples=analysis.samples[:2])
        assert fs.sum() < analysis.matrix_raw.shape[0]


def test_measure_coverage(various_analysis):
    for analysis in various_analysis:
        with pytest.raises(IOError):
            analysis.measure_coverage()


def test_consensus_set_loading(various_analysis):
    for analysis in various_analysis:
        assert hasattr(analysis, "sites")
        assert isinstance(analysis.sites, pybedtools.BedTool)


def test_coverage_matrix_loading(various_analysis):
    for analysis in various_analysis:
        assert hasattr(analysis, "matrix_raw")
        assert isinstance(analysis.matrix_raw, pd.DataFrame)
        assert analysis.matrix_raw.dtypes.all() == int


def test_set_consensus_set(various_analysis):
    for analysis in various_analysis:
        peaks = os.path.join(
            analysis.results_dir, analysis.name + ".peak_set.bed")
        analysis.set_consensus_sites(peaks)
        assert hasattr(analysis, "sites")
        sites = pd.read_csv(peaks, header=None)
        assert len(analysis.sites) == sites.shape[0]


def test_rpm_normalization(various_analysis):
    for analysis in various_analysis:
        qnorm = analysis.normalize_rpm(save=False)
        assert hasattr(analysis, "matrix_norm")
        assert isinstance(qnorm, pd.DataFrame)
        assert qnorm.dtypes.all() == np.float
        assert qnorm.isnull().sum().sum() == 0
        rpm_file = os.path.join(
            analysis.results_dir, analysis.name + ".matrix_norm.csv"
        )
        assert not file_exists(rpm_file)
        qnorm = analysis.normalize_rpm(save=True)
        assert hasattr(analysis, "matrix_norm")
        assert isinstance(qnorm, pd.DataFrame)
        assert qnorm.dtypes.all() == np.float
        assert qnorm.isnull().sum().sum() == 0
        assert file_exists_and_not_empty(rpm_file)
        assert hasattr(analysis, "matrix_norm")


def test_python_quantile_normalization(various_analysis):
    for analysis in various_analysis:
        f = os.path.join(
            analysis.results_dir, analysis.name + ".matrix_norm.csv")
        qnorm_p = analysis.normalize_quantiles(implementation="Python", save=1)
        assert hasattr(analysis, "matrix_norm")
        assert isinstance(qnorm_p, pd.DataFrame)
        assert qnorm_p.dtypes.all() == np.float
        assert qnorm_p.isnull().sum().sum() == 0
        assert file_exists_and_not_empty(f)


@pytest.mark.skipif(
    not R,
    reason=R_REASON)
def test_r_quantile_normalization(various_analysis):
    for analysis in various_analysis:
        f = os.path.join(
            analysis.results_dir, analysis.name + ".matrix_norm.csv")
        qnorm_r = analysis.normalize_quantiles(implementation="R", save=True)
        assert hasattr(analysis, "matrix_norm")
        assert isinstance(qnorm_r, pd.DataFrame)
        assert qnorm_r.dtypes.all() == np.float
        assert qnorm_r.isnull().sum().sum() == 0
        assert file_exists_and_not_empty(f)


# def test_quantile_normalization(various_analysis):
#     for analysis in various_analysis:
#         f = os.path.join(
#             analysis.results_dir, analysis.name + ".matrix_norm.csv")
#         qnorm_p = analysis.normalize_quantiles(implementation="Python", save=1)
#         assert hasattr(analysis, "matrix_norm")
#         assert isinstance(qnorm_p, pd.DataFrame)
#         assert qnorm_p.dtypes.all() == np.float
#         assert qnorm_p.isnull().sum().sum() == 0
#         assert file_exists_and_not_empty(f)

#         qnorm_r = analysis.normalize_quantiles(implementation="R", save=True)
#         assert hasattr(analysis, "matrix_norm")
#         assert isinstance(qnorm_r, pd.DataFrame)
#         assert qnorm_r.dtypes.all() == np.float
#         assert qnorm_r.isnull().sum().sum() == 0
#         assert file_exists_and_not_empty(f)

#         cors = list()
#         for col in qnorm_p.columns:
#             cors.append(scipy.stats.pearsonr(qnorm_p[col], qnorm_r[col])[0])
#         assert all(np.array(cors) > 0.99)


@pytest.mark.skipif(CI, reason="CQN normalization not testable in CI")
def test_cqn_normalization(atac_analysis):
    qnorm = atac_analysis.normalize_cqn()
    assert hasattr(atac_analysis, "matrix_norm")
    assert isinstance(qnorm, pd.DataFrame)
    assert qnorm.dtypes.all() == np.float
    assert qnorm.isnull().sum().sum() == 0
    file = os.path.join(
        atac_analysis.results_dir, atac_analysis.name + ".matrix_norm.csv")
    assert file_exists_and_not_empty(file)


def test_pca_normalization(atac_analysis):
    qnorm = atac_analysis.normalize_pca(pc=1)
    assert hasattr(atac_analysis, "matrix_norm")
    assert isinstance(qnorm, pd.DataFrame)
    assert qnorm.dtypes.all() == np.float
    assert qnorm.isnull().sum().sum() == 0
    file = os.path.join(
        atac_analysis.results_dir, atac_analysis.name + ".matrix_norm.csv")
    assert file_exists_and_not_empty(file)
    plot = os.path.join(atac_analysis.results_dir, "PCA_based_batch_correction.svg")
    assert file_exists_and_not_empty(plot)


def test_vst_normalization(atac_analysis):
    qnorm = atac_analysis.normalize_vst(fitType="mean")
    assert hasattr(atac_analysis, "matrix_norm")
    assert isinstance(qnorm, pd.DataFrame)
    assert qnorm.dtypes.all() == np.float
    assert qnorm.isnull().sum().sum() == 0
    file = os.path.join(
        atac_analysis.results_dir, atac_analysis.name + ".matrix_norm.csv")
    assert file_exists_and_not_empty(file)


def test_normalize(atac_analysis):
    rpmnorm = atac_analysis.normalize_rpm(save=False)
    assert isinstance(rpmnorm, pd.DataFrame)
    assert hasattr(atac_analysis, "matrix_norm")
    del atac_analysis.matrix_norm

    rpmnorm_d = atac_analysis.normalize(method="rpm", save=False)
    assert isinstance(rpmnorm_d, pd.DataFrame)
    assert hasattr(atac_analysis, "matrix_norm")
    assert np.array_equal(rpmnorm_d, rpmnorm)
    del atac_analysis.matrix_norm

    qnorm = atac_analysis.normalize_quantiles(save=False)
    assert hasattr(atac_analysis, "matrix_norm")
    del atac_analysis.matrix_norm

    qnorm_d = atac_analysis.normalize(method="quantile", save=False)
    assert isinstance(qnorm_d, pd.DataFrame)
    assert hasattr(atac_analysis, "matrix_norm")
    assert np.array_equal(qnorm_d, qnorm)
    del atac_analysis.matrix_norm

    if not CI:
        norm = atac_analysis.normalize(method="cqn", save=False)
        assert isinstance(norm, pd.DataFrame)
        assert hasattr(atac_analysis, "matrix_norm")

    norm = atac_analysis.normalize(method="pca", pc=1, save=False)
    assert hasattr(atac_analysis, "matrix_norm")
    del atac_analysis.matrix_norm

    norm = atac_analysis.normalize(method="vst", save=False, fitType="mean")
    assert hasattr(atac_analysis, "matrix_norm")
    del atac_analysis.matrix_norm


def test_get_matrix_stats(various_analysis):
    for analysis in various_analysis:
        annot = analysis.get_matrix_stats(matrix="matrix_raw")
        output = os.path.join(
            analysis.results_dir, analysis.name + ".stats_per_feature.csv"
        )
        assert file_exists_and_not_empty(output)
        assert isinstance(annot, pd.DataFrame)
        cols = ["mean", "variance", "std_deviation", "dispersion", "qv2", "amplitude"]
        assert all([x in annot.columns.tolist() for x in cols])


def test_get_peak_gene_annotation(atac_analysis):
    reference_dir = ATACSeqAnalysis._format_string_with_environment_variables(
        _CONFIG["preferences"]["root_reference_dir"]
    )
    if reference_dir is None:
        reference_dir = os.path.join(atac_analysis.root_dir, "reference")

    mapping = {"hg19": "grch37", "hg38": "grch38", "mm10": "grcm38"}

    annot = atac_analysis.get_peak_gene_annotation(max_dist=1e10)
    tss = os.path.join(
        reference_dir,
        "{}.{}.gene_annotation.protein_coding.tss.bed".format(
            atac_analysis.organism, mapping[atac_analysis.genome]
        ),
    )
    assert file_exists_and_not_empty(tss)
    assert isinstance(annot, pd.DataFrame)
    assert annot.shape[0] >= len(atac_analysis.sites)


def test_get_peak_genomic_location(atac_analysis):
    reference_dir = ATACSeqAnalysis._format_string_with_environment_variables(
        _CONFIG["preferences"]["root_reference_dir"]
    )
    if reference_dir is None:
        reference_dir = os.path.join(atac_analysis.root_dir, "reference")
    prefix = os.path.join(reference_dir, "{}.{}.genomic_context")

    fs = [
        prefix + a
        for a in [
            ".bed",
            ".exon.bed",
            ".genebody.bed",
            ".intergenic.bed",
            ".intron.bed",
            ".promoter.bed",
            ".utr3.bed",
            ".utr5.bed",
        ]
    ]
    annot = atac_analysis.get_peak_genomic_location()

    # check annotation files are produced
    mapping = {"hg19": "grch37", "hg38": "grch38", "mm10": "grcm38"}
    for f in fs:
        f = f.format(atac_analysis.organism, mapping[atac_analysis.genome])
        assert file_exists_and_not_empty(f)

    assert isinstance(annot, pd.DataFrame)
    assert annot.shape[0] >= len(atac_analysis.sites)


def test_peak_chromatin_state(atac_analysis, chrom_file):
    prefix = os.path.join(atac_analysis.results_dir, atac_analysis.name)
    fs = [
        prefix + ".chrom_state_annotation.csv",
        prefix + ".chrom_state_annotation_mapping.csv",
        prefix + ".chrom_state_annotation_background.csv",
        prefix + ".chrom_state_annotation_background_mapping.csv",
    ]

    attrs = [
        "chrom_state_annotation",
        "chrom_state_annotation_b",
        "chrom_state_annotation_mapping",
        "chrom_state_annotation_b_mapping",
    ]

    annot = atac_analysis.get_peak_chromatin_state(chrom_state_file=chrom_file)
    assert isinstance(annot, pd.DataFrame)
    assert annot.shape[0] >= len(atac_analysis.sites)

    for f in fs:
        assert file_exists_and_not_empty(f)
    for attr in attrs:
        assert hasattr(atac_analysis, attr)


def test_annotate_features(atac_analysis, chrom_file):
    atac_analysis.get_peak_chromatin_state(chrom_state_file=chrom_file)
    atac_analysis.get_matrix_stats(matrix="matrix_raw")
    atac_analysis.get_peak_gene_annotation(max_dist=1e10)
    atac_analysis.get_peak_genomic_location()
    atac_analysis.annotate_features(matrix="matrix_raw")
    f = os.path.join(atac_analysis.results_dir, atac_analysis.name + ".matrix_features.csv")
    assert hasattr(atac_analysis, "matrix_features")
    assert file_exists_and_not_empty(f)

    cols = [
        "gene_name",
        "strand",
        "distance",  # from gene_annotation
        "mean",
        "variance",
        "std_deviation",
        "dispersion",
        "qv2",
        "amplitude",
        "iqr",
    ]  # from stats

    cols += ["genomic_region"]  # from genomic_location

    cols += ["chromatin_state"]  # from chromatin_state

    assert all([c in atac_analysis.matrix_features.columns.tolist() for c in cols])


def test_plot_raw_coverage(various_analysis):
    for analysis in various_analysis:
        analysis.plot_raw_coverage()
        output = os.path.join(
            analysis.results_dir, analysis.name + ".raw_counts.violinplot.svg"
        )
        assert file_exists_and_not_empty(output)

        attr = "A"
        analysis.plot_raw_coverage(by_attribute=attr)
        output = os.path.join(
            analysis.results_dir,
            analysis.name + ".raw_counts.violinplot.by_{}.svg".format(attr),
        )
        assert file_exists_and_not_empty(output)


def test_get_sex_chrom_ratio(analysis_normalized):
    fs = [
        ".sex_chrom_ratio.csv",
        ".sex_chrom_ratio.clustermap.svg",
        ".sex_chrom_ratio.rank_vs_ratio.svg"]
    with analysis_normalized as a:
        fs = [os.path.join(a.results_dir, a.name + f) for f in fs]
        a.get_sex_chrom_ratio()

        for f in fs:
            assert file_exists_and_not_empty(f)


def test_get_sex_chrom_ratio_wrong_sex_chroms(analysis_normalized):
    with analysis_normalized as a:
        with pytest.raises(ValueError):
            a.get_sex_chrom_ratio(sex_chroms=['qwe13213'])


def test_get_sex_chrom_ratio_no_sex_chroms(analysis_normalized):
    with analysis_normalized as a:
        matrix = a.matrix_norm.loc[~a.matrix_norm.index.str.contains("chrX|chrY")]
        with pytest.raises(ValueError):
            a.get_sex_chrom_ratio(matrix)


@pytest.fixture
def peak_outputs(atac_analysis_with_input_files):
    outputs = [
        ".peak_location.per_sample.svg",
        ".lengths.svg",
        ".peak_lengths.per_sample.svg",
        ".total_open_chromatin_space.per_sample.svg",
        ".open_chromatin_space.csv"]
    with atac_analysis_with_input_files as a:
        outputs = [os.path.join(a.results_dir, "peak_characteristics", a.name + f) for f in outputs]

    return outputs


class Test_plot_peak_characteristics:
    # def test_no_arguments(
    #         self, atac_analysis_with_input_files, peak_outputs):
    #     with atac_analysis_with_input_files as a:
    #         a.plot_peak_characteristics()

    #         for f in peak_outputs:
    #             assert file_exists_and_not_empty(f)

    # def test_plot_peak_characteristics_by_group(
    #         self, atac_analysis_with_input_files, peak_outputs):
    #     with atac_analysis_with_input_files as a:
    #         g = a.group_attributes[-1]
    #         outputs = [
    #             ".peak_lengths.per_{}.svg".format(g),
    #             ".total_open_chromatin_space.per_{}.svg".format(g)]
    #         outputs = [os.path.join(a.results_dir, "peak_characteristics", a.name + f) for f in outputs]

    #         a.plot_peak_characteristics(by_attribute=g)

    #         for f in peak_outputs + outputs:
    #             assert file_exists_and_not_empty(f)

    # def test_plot_peak_characteristics_with_closest_tss(
    #         self, atac_analysis_with_input_files, peak_outputs):
    #     with atac_analysis_with_input_files as a:
    #         outputs = [
    #             ".tss_distance.svg"]
    #         outputs = [os.path.join(a.results_dir, "peak_characteristics", a.name + f) for f in outputs]

    #         a.get_peak_gene_annotation()
    #         a.plot_peak_characteristics()

    #         for f in peak_outputs + outputs:
    #             assert file_exists_and_not_empty(f)

    # def test_plot_peak_characteristics_with_genomic_location(
    #         self, atac_analysis_with_input_files, peak_outputs):
    #     with atac_analysis_with_input_files as a:
    #         outputs = [
    #             ".genomic_regions.svg"]
    #         outputs = [os.path.join(a.results_dir, "peak_characteristics", a.name + f) for f in outputs]

    #         a.get_peak_genomic_location()
    #         a.plot_peak_characteristics()

    #         for f in peak_outputs + outputs:
    #             assert file_exists_and_not_empty(f)

    # def test_plot_peak_characteristics_with_chromatin_state(
    #         self, atac_analysis_with_input_files, chrom_file, peak_outputs):
    #     with atac_analysis_with_input_files as a:
    #         outputs = [
    #             ".chromatin_states.svg"]
    #         outputs = [os.path.join(a.results_dir, "peak_characteristics", a.name + f) for f in outputs]

    #         a.get_peak_chromatin_state(chrom_file)
    #         a.plot_peak_characteristics()

    #         for f in peak_outputs + outputs:
    #             assert file_exists_and_not_empty(f)

    # def test_plot_peak_characteristics_with_both(
    #         self, atac_analysis_with_input_files, chrom_file, peak_outputs):
    #     with atac_analysis_with_input_files as a:
    #         outputs = [
    #             ".genomic_regions.svg",
    #             ".chromatin_states.svg",
    #             ".genomic_region_and_chromatin_states.svg"]
    #         outputs = [os.path.join(a.results_dir, "peak_characteristics", a.name + f) for f in outputs]

    #         a.get_peak_genomic_location()
    #         a.get_peak_chromatin_state(chrom_file)
    #         a.plot_peak_characteristics()

    #         for f in peak_outputs + outputs:
    #             assert file_exists_and_not_empty(f)

    # def test_plot_peak_characteristics_with_both_and_stats(
    #         self, atac_analysis_with_input_files, chrom_file, peak_outputs):
    #     with atac_analysis_with_input_files as a:
    #         outputs = [
    #             ".genomic_regions.svg",
    #             ".chromatin_states.svg",
    #             ".genomic_region_and_chromatin_states.svg",
    #             ".mean_vs_iqr.svg",
    #             ".mean_vs_amplitude.svg",
    #             ".mean_vs_qv2.svg",
    #             ".mean_vs_dispersion.svg",
    #             ".mean_vs_variance.svg",
    #             ".mean_vs_std_deviation.svg",
    #             ".qv2.distplot.svg",
    #             ".iqr.distplot.svg",
    #             ".amplitude.distplot.svg",
    #             ".variance.distplot.svg",
    #             ".std_deviation.distplot.svg",
    #             ".mean.distplot.svg",
    #             ".dispersion.distplot.svg"]
    #         outputs = [os.path.join(a.results_dir, "peak_characteristics", a.name + f) for f in outputs]

    #         a.get_peak_genomic_location()
    #         a.get_peak_chromatin_state(chrom_file)
    #         a.get_matrix_stats()
    #         a.plot_peak_characteristics()

    #         for f in peak_outputs + outputs:
    #             assert file_exists_and_not_empty(f)

    def test_plot_peak_characteristics_with_both_and_stats_and_support(
            self, atac_analysis_with_input_files, chrom_file, peak_outputs):
        with atac_analysis_with_input_files as a:
            outputs = [
                ".genomic_regions.svg",
                ".chromatin_states.svg",
                ".genomic_region_and_chromatin_states.svg",
                ".mean_vs_iqr.svg",
                ".mean_vs_amplitude.svg",
                ".mean_vs_qv2.svg",
                ".mean_vs_dispersion.svg",
                ".mean_vs_variance.svg",
                ".mean_vs_std_deviation.svg",
                ".qv2.distplot.svg",
                ".iqr.distplot.svg",
                ".amplitude.distplot.svg",
                ".variance.distplot.svg",
                ".std_deviation.distplot.svg",
                ".mean.distplot.svg",
                ".dispersion.distplot.svg",
                ".mean_vs_support.svg"]
            outputs = [os.path.join(a.results_dir, "peak_characteristics", a.name + f) for f in outputs]

            a.calculate_peak_support()
            a.get_peak_genomic_location()
            a.get_peak_chromatin_state(chrom_file)
            a.get_matrix_stats()
            a.plot_peak_characteristics()

            for f in peak_outputs + outputs:
                assert file_exists_and_not_empty(f)
