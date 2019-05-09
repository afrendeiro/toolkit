#!/usr/bin/env python


import os

from .data_generator import generate_project
from ngs_toolkit import _CONFIG
from ngs_toolkit.atacseq import ATACSeqAnalysis
from ngs_toolkit.utils import download_gzip_file
import numpy as np
import pandas as pd
from peppy import Project
import pybedtools
import pytest


travis = "TRAVIS" in os.environ


@pytest.fixture
def various_analysis(tmp_path):
    tmp_path = str(tmp_path)  # for Python2

    # Let's make several "reallish" test projects
    to_test = list()
    project_prefix_name = "test-project"
    data_type = "ATAC-seq"
    genome_assemblies = [("human", "hg19"), ("human", "hg38"), ("mouse", "mm10")]

    n_factors = [1, 2, 3][0]
    n_variables = [100, 1000, 10000][0]
    n_replicates = [1, 2, 5][0]
    for organism, genome_assembly in genome_assemblies:
        project_name = "{}_{}_{}_{}_{}_{}".format(
            project_prefix_name,
            data_type,
            genome_assembly,
            n_factors,
            n_variables,
            n_replicates,
        )

        generate_project(
            output_dir=tmp_path,
            project_name=project_name,
            organism=organism,
            genome_assembly=genome_assembly,
            data_type=data_type,
            n_factors=n_factors,
            n_replicates=n_replicates,
            n_variables=n_variables,
        )

        # first edit the defaul path to the annotation sheet
        config = os.path.join(tmp_path, project_name, "metadata", "project_config.yaml")
        prj_path = os.path.join(tmp_path, project_name)
        os.chdir(prj_path)

        # project and associated analysis
        analysis = ATACSeqAnalysis(
            name=project_name,
            prj=Project(config),
            results_dir=os.path.join(prj_path, "results"),
        )
        analysis.load_data()

        to_test.append(analysis)
    return to_test


@pytest.fixture
def analysis(tmp_path):
    tmp_path = str(tmp_path)  # for Python2

    # Let's make several "reallish" test projects
    project_prefix_name = "test-project"
    data_type = "ATAC-seq"
    organism, genome_assembly = ("human", "hg38")

    n_factors = [1, 2, 3][0]
    n_variables = [100, 1000, 10000][0]
    n_replicates = [1, 2, 5][0]
    project_name = "{}_{}_{}_{}_{}_{}".format(
        project_prefix_name,
        data_type,
        genome_assembly,
        n_factors,
        n_variables,
        n_replicates,
    )

    generate_project(
        output_dir=tmp_path,
        project_name=project_name,
        organism=organism,
        genome_assembly=genome_assembly,
        data_type=data_type,
        n_factors=n_factors,
        n_replicates=n_replicates,
        n_variables=n_variables,
    )

    # first edit the defaul path to the annotation sheet
    config = os.path.join(tmp_path, project_name, "metadata", "project_config.yaml")
    prj_path = os.path.join(tmp_path, project_name)
    os.chdir(prj_path)

    # project and associated analysis
    analysis = ATACSeqAnalysis(
        name=project_name,
        prj=Project(config),
        results_dir=os.path.join(prj_path, "results"),
    )
    analysis.load_data()

    return analysis


@pytest.fixture
def chrom_file():
    url = (
        "https://egg2.wustl.edu/roadmap/data/byFileType/"
        + "chromhmmSegmentations/ChmmModels/coreMarks/jointModel/"
        + "final/E002_15_coreMarks_hg38lift_dense.bed.gz"
    )
    chrom_state_file = os.path.abspath("E002_15_coreMarks_hg38lift_dense.bed")
    download_gzip_file(url, chrom_state_file)
    return chrom_state_file


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
def test_cqn_normalization(analysis):
    # At some point, downloading a genome reference in Travis
    # caused memory error.
    # This should now be fixed by implementing download/decompressing
    # functions working in chunks
    try:
        qnorm = analysis.normalize_cqn()
    except OSError:
        if travis:
            pytest.skip()
        else:
            raise
    assert qnorm.dtypes.all() == np.float
    file = os.path.join(analysis.results_dir, analysis.name + ".matrix_norm.csv")
    assert os.path.exists(file)
    assert os.stat(file).st_size > 0


def test_normalize(analysis):
    qnorm = analysis.normalize_rpm(save=False)
    assert isinstance(qnorm, pd.DataFrame)
    assert hasattr(analysis, "matrix_norm")
    del analysis.matrix_norm

    qnorm_d = analysis.normalize(method="rpm", save=False)
    assert isinstance(qnorm_d, pd.DataFrame)
    assert hasattr(analysis, "matrix_norm")
    assert np.array_equal(qnorm_d, qnorm)
    del analysis.matrix_norm

    qnorm = analysis.normalize_quantiles(save=False)
    assert hasattr(analysis, "matrix_norm")
    del analysis.matrix_norm

    qnorm_d = analysis.normalize(method="quantile", save=False)
    assert isinstance(qnorm_d, pd.DataFrame)
    assert hasattr(analysis, "matrix_norm")
    assert np.array_equal(qnorm_d, qnorm)
    del analysis.matrix_norm

    # At some point, downloading a genome reference in Travis
    # caused memory error.
    # This should now be fixed by implementing download/decompressing
    # functions working in chunks
    if not travis:
        norm = analysis.normalize(method="cqn", save=False)
        assert isinstance(norm, pd.DataFrame)
        assert hasattr(analysis, "matrix_norm")


def test_get_matrix_stats(various_analysis):
    for analysis in various_analysis:
        annot = analysis.get_matrix_stats(matrix="matrix_raw")
        output = os.path.join(
            analysis.results_dir, analysis.name + ".stats_per_region.csv"
        )
        assert os.path.exists(output)
        assert os.stat(output).st_size > 0
        assert isinstance(annot, pd.DataFrame)
        cols = ["mean", "variance", "std_deviation", "dispersion", "qv2", "amplitude"]
        assert all([x in annot.columns.tolist() for x in cols])


def test_get_peak_gene_annotation(analysis):
    reference_dir = ATACSeqAnalysis._format_string_with_environment_variables(
        _CONFIG["preferences"]["root_reference_dir"]
    )
    if reference_dir is None:
        reference_dir = os.path.join(analysis.root_dir, "reference")

    mapping = {"hg19": "grch37", "hg38": "grch38", "mm10": "grcm38"}

    annot = analysis.get_peak_gene_annotation(max_dist=1e10)
    tss = os.path.join(
        reference_dir,
        "{}.{}.gene_annotation.protein_coding.tss.bed".format(
            analysis.organism, mapping[analysis.genome]
        ),
    )
    assert os.path.exists(tss)
    assert os.stat(tss).st_size > 0
    assert isinstance(annot, pd.DataFrame)
    assert annot.shape[0] >= len(analysis.sites)


def test_get_peak_genomic_location(analysis):
    reference_dir = ATACSeqAnalysis._format_string_with_environment_variables(
        _CONFIG["preferences"]["root_reference_dir"]
    )
    if reference_dir is None:
        reference_dir = os.path.join(analysis.root_dir, "reference")
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
        annot = analysis.get_peak_genomic_location()
    except OSError:
        if travis:
            pytest.skip()
        else:
            raise

    # check annotation files are produced
    mapping = {"hg19": "grch37", "hg38": "grch38", "mm10": "grcm38"}
    for f in fs:
        f = f.format(analysis.organism, mapping[analysis.genome])
        assert os.path.exists(f)
        assert os.stat(f).st_size > 0

    assert isinstance(annot, pd.DataFrame)
    assert annot.shape[0] >= len(analysis.sites)


def test_peak_chromatin_state(analysis, chrom_file):
    prefix = os.path.join(analysis.results_dir, analysis.name)
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

    annot = analysis.get_peak_chromatin_state(chrom_state_file=chrom_file)
    assert isinstance(annot, pd.DataFrame)
    assert annot.shape[0] >= len(analysis.sites)

    for f in fs:
        assert os.path.exists(f)
        assert os.stat(f).st_size > 0
    for attr in attrs:
        assert hasattr(analysis, attr)


def test_annotate(analysis, chrom_file):
    analysis.get_peak_chromatin_state(chrom_state_file=chrom_file)
    analysis.get_matrix_stats(matrix="matrix_raw")
    analysis.get_peak_gene_annotation(max_dist=1e10)
    # At some point, downloading a genome reference in Travis
    # caused memory error.
    # This should now be fixed by implementing download/decompressing
    # functions working in chunks
    failed = False
    try:
        analysis.get_peak_genomic_location()
    except OSError:
        if travis:
            failed = True
        else:
            raise
    analysis.annotate(matrix="matrix_raw")
    f = os.path.join(analysis.results_dir, analysis.name + ".matrix_features.csv")
    assert hasattr(analysis, "matrix_features")
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

    assert all([c in analysis.matrix_features.columns.tolist() for c in cols])


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
