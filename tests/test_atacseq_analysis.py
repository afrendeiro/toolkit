#!/usr/bin/env python


import os

from ngs_toolkit import _CONFIG
from ngs_toolkit.atacseq import ATACSeqAnalysis
import numpy as np
import pandas as pd
import pybedtools
import pytest


travis = "TRAVIS" in os.environ


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
        peaks = os.path.join(analysis.results_dir, analysis.name + ".peak_set.bed")
        analysis.set_consensus_sites(peaks)
        assert hasattr(analysis, "sites")
        sites = pd.read_csv(peaks, header=None)
        assert len(analysis.sites) == sites.shape[0]


def test_rpm_normalization(various_analysis):
    for analysis in various_analysis:
        qnorm = analysis.normalize_rpm(save=False)
        assert qnorm.dtypes.all() == np.float
        assert hasattr(analysis, "matrix_norm")
        rpm_file = os.path.join(
            analysis.results_dir, analysis.name + ".matrix_norm.csv"
        )
        assert not os.path.exists(rpm_file)
        qnorm = analysis.normalize_rpm(save=True)
        assert os.path.exists(rpm_file)
        assert os.stat(rpm_file).st_size > 0
        assert hasattr(analysis, "matrix_norm")


def test_quantile_normalization(various_analysis):
    for analysis in various_analysis:
        f = os.path.join(analysis.results_dir, analysis.name + ".matrix_norm.csv")
        qnorm_p = analysis.normalize_quantiles(implementation="Python", save=True)
        assert isinstance(qnorm_p, pd.DataFrame)
        assert hasattr(analysis, "matrix_norm")
        assert os.path.exists(f)
        assert os.stat(f).st_size > 0
        del analysis.matrix_norm
        os.remove(f)

        qnorm_r = analysis.normalize_quantiles(implementation="R", save=True)
        assert isinstance(qnorm_r, pd.DataFrame)
        assert hasattr(analysis, "matrix_norm")
        assert os.path.exists(f)
        assert os.stat(f).st_size > 0

        # cors = list()
        # for col in qnorm_p.columns:
        #     cors.append(scipy.stats.pearsonr(qnorm_p[col], qnorm_r[col])[0])
        # assert all(np.array(cors) > 0.99)


@pytest.mark.skipif(travis, reason="This is anyway tested after")
def test_cqn_normalization(atac_analysis):
    # At some point, downloading a genome reference in Travis
    # caused memory error.
    # This should now be fixed by implementing download/decompressing
    # functions working in chunks
    try:
        qnorm = atac_analysis.normalize_cqn()
    except OSError:
        if travis:
            pytest.skip()
        else:
            raise
    assert qnorm.dtypes.all() == np.float
    file = os.path.join(atac_analysis.results_dir, atac_analysis.name + ".matrix_norm.csv")
    assert os.path.exists(file)
    assert os.stat(file).st_size > 0


def test_normalize(atac_analysis):
    qnorm = atac_analysis.normalize_rpm(save=False)
    assert isinstance(qnorm, pd.DataFrame)
    assert hasattr(atac_analysis, "matrix_norm")
    del atac_analysis.matrix_norm

    qnorm_d = atac_analysis.normalize(method="rpm", save=False)
    assert isinstance(qnorm_d, pd.DataFrame)
    assert hasattr(atac_analysis, "matrix_norm")
    assert np.array_equal(qnorm_d, qnorm)
    del atac_analysis.matrix_norm

    qnorm = atac_analysis.normalize_quantiles(save=False)
    assert hasattr(atac_analysis, "matrix_norm")
    del atac_analysis.matrix_norm

    qnorm_d = atac_analysis.normalize(method="quantile", save=False)
    assert isinstance(qnorm_d, pd.DataFrame)
    assert hasattr(atac_analysis, "matrix_norm")
    assert np.array_equal(qnorm_d, qnorm)
    del atac_analysis.matrix_norm

    # At some point, downloading a genome reference in Travis
    # caused memory error.
    # This should now be fixed by implementing download/decompressing
    # functions working in chunks
    if not travis:
        norm = atac_analysis.normalize(method="cqn", save=False)
        assert isinstance(norm, pd.DataFrame)
        assert hasattr(atac_analysis, "matrix_norm")


def test_get_matrix_stats(various_analysis):
    for analysis in various_analysis:
        annot = analysis.get_matrix_stats(matrix="matrix_raw")
        output = os.path.join(
            analysis.results_dir, analysis.name + ".stats_per_feature.csv"
        )
        assert os.path.exists(output)
        assert os.stat(output).st_size > 0
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
    assert os.path.exists(tss)
    assert os.stat(tss).st_size > 0
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
    # At some point, downloading a genome reference in Travis
    # caused memory error.
    # This should now be fixed by implementing download/decompressing
    # functions working in chunks
    try:
        annot = atac_analysis.get_peak_genomic_location()
    except OSError:
        if travis:
            pytest.skip()
        else:
            raise

    # check annotation files are produced
    mapping = {"hg19": "grch37", "hg38": "grch38", "mm10": "grcm38"}
    for f in fs:
        f = f.format(atac_analysis.organism, mapping[atac_analysis.genome])
        assert os.path.exists(f)
        assert os.stat(f).st_size > 0

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
        assert os.path.exists(f)
        assert os.stat(f).st_size > 0
    for attr in attrs:
        assert hasattr(atac_analysis, attr)


def test_annotate_features(atac_analysis, chrom_file):
    atac_analysis.get_peak_chromatin_state(chrom_state_file=chrom_file)
    atac_analysis.get_matrix_stats(matrix="matrix_raw")
    atac_analysis.get_peak_gene_annotation(max_dist=1e10)
    # At some point, downloading a genome reference in Travis
    # caused memory error.
    # This should now be fixed by implementing download/decompressing
    # functions working in chunks
    failed = False
    try:
        atac_analysis.get_peak_genomic_location()
    except OSError:
        if travis:
            failed = True
        else:
            raise
    atac_analysis.annotate_features(matrix="matrix_raw")
    f = os.path.join(atac_analysis.results_dir, atac_analysis.name + ".matrix_features.csv")
    assert hasattr(atac_analysis, "matrix_features")
    assert os.path.exists(f)
    assert os.stat(f).st_size > 0

    cols = [
        "gene_name",
        "strand",
        "distance",  # from gene_annotation
        "chromatin_state",  # from chromatin_state
        "mean",
        "variance",
        "std_deviation",
        "dispersion",
        "qv2",
        "amplitude",
        "iqr",
    ]  # from stats

    if not failed:
        cols += ["genomic_region"]  # from genomic_location

    assert all([c in atac_analysis.matrix_features.columns.tolist() for c in cols])


def test_plot_raw_coverage(various_analysis):
    for analysis in various_analysis:
        analysis.plot_raw_coverage()
        output = os.path.join(
            analysis.results_dir, analysis.name + ".raw_counts.violinplot.svg"
        )
        assert os.path.exists(output)
        assert os.stat(output).st_size > 0

        attr = "a"
        analysis.plot_raw_coverage(by_attribute=attr)
        output = os.path.join(
            analysis.results_dir,
            analysis.name + ".raw_counts.violinplot.by_{}.svg".format(attr),
        )
        assert os.path.exists(output)
        assert os.stat(output).st_size > 0
