#!/usr/bin/env python


import os

from .data_generator import generate_project
from ngs_toolkit import Analysis, ATACSeqAnalysis
from peppy import Project
import pytest


travis = "TRAVIS" in os.environ


# TODO: test having no config set
# TODO: test differential analysis with many factors
# TODO: test subproject initialization


@pytest.fixture
def empty_analysis():
    return ATACSeqAnalysis()


@pytest.fixture
def null_analysis():
    analysis = ATACSeqAnalysis()
    analysis.organism = None
    analysis.genome = None
    analysis.sites = None
    return analysis


@pytest.fixture
def full_analysis():
    analysis = ATACSeqAnalysis()
    analysis.organism = "human"
    analysis.genome = "hg38"
    analysis.sites = "hg38"
    analysis.samples = []
    return analysis


@pytest.fixture
def analysis(tmp_path):
    return Analysis(name="test-project")


@pytest.fixture
def atac_analysis(tmp_path):
    tmp_path = str(tmp_path)  # for Python2

    # Let's make several "reallish" test projects
    project_prefix_name = "test-project"
    data_type = "ATAC-seq"
    organism, genome_assembly = ("human", "hg38")

    n_factors = 1
    n_variables = 1000
    n_replicates = 3
    project_name = "_".join(
        str(x) for x in [
            project_prefix_name,
            data_type,
            genome_assembly,
            n_factors,
            n_variables,
            n_replicates])

    generate_project(
        output_dir=tmp_path,
        project_name=project_name,
        organism=organism,
        genome_assembly=genome_assembly,
        data_type=data_type,
        n_factors=n_factors,
        n_replicates=n_replicates,
        n_variables=n_variables)

    # first edit the defaul path to the annotation sheet
    config = os.path.join(
        tmp_path, project_name, "metadata", "project_config.yaml")
    prj_path = os.path.join(
        tmp_path, project_name)

    # project and associated analysis
    analysis = ATACSeqAnalysis(
        name=project_name,
        prj=Project(config),
        results_dir=os.path.join(prj_path, "results"))
    analysis.load_data()

    return analysis


@pytest.fixture
def analysis_normalized(atac_analysis):
    atac_analysis.normalize(method="rpm")
    return atac_analysis


@pytest.fixture
def analysis_annotated(analysis_normalized):
    analysis_normalized.annotate_samples()
    return analysis_normalized


@pytest.fixture
def analysis_with_differential(analysis_normalized):
    analysis_normalized.get_peak_gene_annotation()
    analysis_normalized.annotate_features()
    analysis_normalized.annotate_samples()
    analysis_normalized.differential_analysis(filter_support=False)
    return analysis_normalized


@pytest.fixture
def analysis_with_differential_enrichment(analysis_with_differential):
    from ngs_toolkit import _CONFIG

    _CONFIG["resources"]["enrichr"]["gene_set_libraries"] = [
        "GO_Biological_Process_2015",
        "NCI-Nature_2016",
    ]
    analysis_with_differential.differential_enrichment(steps=["enrichr"])

    return analysis_with_differential


@pytest.fixture
def atac_analysis_many_factors(tmp_path):
    tmp_path = str(tmp_path)  # for Python2

    # Let's make several "reallish" test projects
    project_prefix_name = "test-project"
    data_type = "ATAC-seq"
    organism, genome_assembly = ("human", "hg38")

    n_factors = 4
    n_variables = 1000
    n_replicates = 4
    project_name = "_".join(
        str(x) for x in [
            project_prefix_name,
            data_type,
            genome_assembly,
            n_factors,
            n_variables,
            n_replicates])

    generate_project(
        output_dir=tmp_path,
        project_name=project_name,
        organism=organism,
        genome_assembly=genome_assembly,
        data_type=data_type,
        n_factors=n_factors,
        n_replicates=n_replicates,
        n_variables=n_variables)

    # first edit the defaul path to the annotation sheet
    config = os.path.join(
        tmp_path, project_name, "metadata", "project_config.yaml")
    prj_path = os.path.join(
        tmp_path, project_name)

    # project and associated analysis
    analysis = ATACSeqAnalysis(
        name=project_name,
        prj=Project(config),
        results_dir=os.path.join(prj_path, "results"))
    analysis.load_data()
    analysis.normalize(method="rpm")
    analysis.annotate_samples()

    return analysis


@pytest.fixture
def rnaseq_analysis(tmp_path):
    tmp_path = str(tmp_path)  # for Python2

    # Let's make several "reallish" test projects
    project_prefix_name = "test-project"
    data_type = "RNA-seq"
    organism, genome_assembly = ("human", "hg38")

    n_factors = 1
    n_variables = 100
    n_replicates = 1
    project_name = "_".join(
        str(x) for x in [
            project_prefix_name,
            data_type,
            genome_assembly,
            n_factors,
            n_variables,
            n_replicates])

    generate_project(
        output_dir=tmp_path,
        project_name=project_name,
        organism=organism,
        genome_assembly=genome_assembly,
        data_type=data_type,
        n_factors=n_factors,
        n_replicates=n_replicates,
        n_variables=n_variables)

    # first edit the defaul path to the annotation sheet
    config = os.path.join(
        tmp_path, project_name, "metadata", "project_config.yaml")
    prj_path = os.path.join(
        tmp_path, project_name)

    # project and associated analysis
    analysis = ATACSeqAnalysis(
        name=project_name,
        prj=Project(config),
        results_dir=os.path.join(prj_path, "results"),
    )
    analysis.load_data()

    return analysis


@pytest.fixture
def various_analysis(tmp_path):
    tmp_path = str(tmp_path)  # for Python2

    # Let's make several "reallish" test projects
    to_test = list()
    project_prefix_name = "test-project"
    data_type = "ATAC-seq"
    genome_assemblies = [("human", "hg19"), ("human", "hg38"), ("mouse", "mm10")]
    factors = [1, 2, 3]
    variables = [100, 1000]  # 10000
    replicates = [1, 2, 5]

    for organism, genome_assembly in genome_assemblies:
        for n_factors in factors:
            for n_variables in variables:
                for n_replicates in replicates:
                    project_name = "_".join(
                        str(x) for x in [
                            project_prefix_name,
                            data_type,
                            genome_assembly,
                            n_factors,
                            n_variables,
                            n_replicates])

                    generate_project(
                        output_dir=tmp_path,
                        project_name=project_name,
                        organism=organism,
                        genome_assembly=genome_assembly,
                        data_type=data_type,
                        n_factors=n_factors,
                        n_replicates=n_replicates,
                        n_variables=n_variables)

                    # first edit the defaul path to the annotation sheet
                    config = os.path.join(
                        tmp_path, project_name, "metadata", "project_config.yaml")
                    prj_path = os.path.join(
                        tmp_path, project_name)

                    # project and associated analysis
                    analysis = ATACSeqAnalysis(
                        name=project_name,
                        prj=Project(config),
                        results_dir=os.path.join(prj_path, "results"),
                    )
                    # make sure the object is loaded with its dataframes
                    analysis.load_data()

                    to_test.append(analysis)
    return to_test


@pytest.fixture
def chrom_file():
    from ngs_toolkit.utils import download_gzip_file
    import pandas as pd

    url = (
        "https://egg2.wustl.edu/roadmap/data/byFileType/"
        + "chromhmmSegmentations/ChmmModels/coreMarks/jointModel/"
        + "final/E002_15_coreMarks_hg38lift_dense.bed.gz"
    )
    chrom_state_file = os.path.abspath("E002_15_coreMarks_hg38lift_dense.bed")
    download_gzip_file(url, chrom_state_file)

    # Test
    assert os.path.exists(chrom_state_file)
    assert os.stat(chrom_state_file).st_size > 0
    b = pd.read_csv(chrom_state_file, skiprows=1, sep="\t")
    assert b.shape == (281837, 9)
    return chrom_state_file
