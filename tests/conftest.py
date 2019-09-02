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


def is_internet_connected(hostname="www.google.com"):
    import socket
    try:
        # see if we can resolve the host name -- tells us if there is
        # a DNS listening
        host = socket.gethostbyname(hostname)
        # connect to the host -- tells us if the host is actually
        # reachable
        s = socket.create_connection((host, 80), 2)
        s.close()
        return True
    except OSError:
        pass
    return False


def file_exists(file):
    import os
    from ngs_toolkit.utils import get_this_file_or_timestamped

    return os.path.exists(get_this_file_or_timestamped(file))


def file_not_empty(file):
    import os
    from ngs_toolkit.utils import get_this_file_or_timestamped

    return os.stat(get_this_file_or_timestamped(file)).st_size > 0


def file_exists_and_not_empty(file):
    import os
    from ngs_toolkit.utils import get_this_file_or_timestamped

    f = get_this_file_or_timestamped(file)

    return os.path.exists(f) and (
        os.stat(f).st_size > 0)


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
def atac_analysis_with_input_files(atac_analysis):
    from ngs_toolkit import _CONFIG
    from .data_generator import generate_sample_input_files

    c = {
        "sample_input_files": {
            "ATAC-seq": {
                "aligned_filtered_bam": "{data_dir}/{sample_name}/mapped/{sample_name}.trimmed.bowtie2.filtered.bam",
                "peaks": "{data_dir}/{sample_name}/peaks/{sample_name}_peaks.narrowPeak",
                "summits": "{data_dir}/{sample_name}/peaks/{sample_name}_summits.bed"},
            "ChIP-seq": {
                "aligned_filtered_bam": "{data_dir}/{sample_name}/mapped/{sample_name}.trimmed.bowtie2.filtered.bam"},
            "CNV": {
                "aligned_filtered_bam": "{data_dir}/{sample_name}/mapped/{sample_name}.trimmed.bowtie2.filtered.bam"},
            "RNA-seq": {
                "aligned_filtered_bam": "{data_dir}/{sample_name}/mapped/{sample_name}.trimmed.bowtie2.filtered.bam",
                "bitseq_counts": "{data_dir}/{sample_name}/bowtie1_{genome}/bitSeq/{sample_name}.counts"}}}
    _CONFIG.update(c)
    atac_analysis.set_samples_input_files()

    generate_sample_input_files(atac_analysis)

    return atac_analysis


@pytest.fixture
def subproject_config(atac_analysis):
    import yaml

    annot = os.path.join(atac_analysis.root_dir, "metadata", "annotation.csv")
    subannot = os.path.join(atac_analysis.root_dir, "metadata", "sample_subannotation.csv")

    yaml_file = os.path.join(atac_analysis.root_dir, "metadata", "project_config.yaml")
    conf = yaml.safe_load(open(yaml_file, 'r'))
    conf['subprojects'] = {"test_subproject": {"metadata": {
        "sample_annotation": annot,
        "sample_subannotation": subannot}}}
    del conf['metadata']['sample_annotation']
    del conf['metadata']['sample_subannotation']

    yaml.safe_dump(conf, open(yaml_file, 'w'))

    return yaml_file


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
    genome_assemblies = [("human", "hg38"), ("mouse", "mm10")]  # ("human", "hg19")
    factors = [1]  # 2, 5
    variables = [100]  # 1000, 10000
    replicates = [1]  # 3, 5

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
