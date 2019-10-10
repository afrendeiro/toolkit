#!/usr/bin/env python


import os

import pytest

from ngs_toolkit import Analysis, ATACSeqAnalysis
from ngs_toolkit.demo import generate_project


CI = ("TRAVIS" in os.environ) or ("GITHUB_WORKFLOW" in os.environ)

CI_NAME = None
BUILD_DIR = os.path.abspath(os.path.curdir)
if CI:
    if ("TRAVIS" in os.environ):
        CI_NAME = "TRAVIS"
        BUILD_DIR = os.environ['TRAVIS_BUILD_DIR']
    elif ("GITHUB_WORKFLOW" in os.environ):
        CI_NAME = "GITHUB"
        BUILD_DIR = os.path.join(
            "home", "runner", "work", "toolkit", "toolkit")

try:
    DEV = os.environ['TRAVIS_BRANCH'] == 'dev'
except KeyError:
    pass
try:
    DEV = os.environ['GITHUB_REF'] == 'dev'
except KeyError:
    import subprocess
    o = subprocess.check_output("git status".split(" "))
    DEV = "dev" in o.decode().split("\n")[0]

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
    tmp_path = str(tmp_path)
    kwargs = {
        'data_type': "ATAC-seq",
        'organism': "human",
        'genome_assembly': "hg38",
        'n_factors': 1,
        'n_features': 250,
        'n_replicates': 2}
    kwargs.update({
        "project_name": "test-project_" + "_".join(
            str(x) for x in kwargs.values()),
        "output_dir": tmp_path})

    return generate_project(**kwargs)


@pytest.fixture
def atac_analysis_with_input_files(tmp_path):
    from ngs_toolkit import _CONFIG
    from ngs_toolkit.demo.data_generator import generate_sample_input_files

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

    tmp_path = str(tmp_path)
    kwargs = {
        'data_type': "ATAC-seq",
        'organism': "human",
        'genome_assembly': "hg38",
        'n_factors': 1,
        'n_features': 10,
        'n_replicates': 2}
    kwargs.update({
        "project_name": "test-project_" + "_".join(
            str(x) for x in kwargs.values()),
        "output_dir": tmp_path})
    an = generate_project(**kwargs)
    an.set_samples_input_files()
    generate_sample_input_files(an, an.matrix_raw)

    return an


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
    analysis_normalized.differential_analysis(
        filter_support=False)
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
    tmp_path = str(tmp_path)
    kwargs = {
        'data_type': "ATAC-seq",
        'organism': "human",
        'genome_assembly': "hg38",
        'n_factors': 4,
        'n_features': 100,
        'n_replicates': 4}
    kwargs.update({
        "project_name": "test-project_" + "_".join(
            str(x) for x in kwargs.values()),
        "output_dir": tmp_path})

    an = generate_project(**kwargs)
    an.load_data()
    an.normalize(method="rpm")
    an.annotate_samples()

    return an


@pytest.fixture
def rnaseq_analysis(tmp_path):
    tmp_path = str(tmp_path)
    kwargs = {
        'data_type': "RNA-seq",
        'organism': "human",
        'genome_assembly': "hg38",
        'n_factors': 1,
        'n_features': 100,
        'n_replicates': 3}
    kwargs.update({
        "project_name": "test-project_" + "_".join(
            str(x) for x in kwargs.values()),
        "output_dir": tmp_path})

    return generate_project(**kwargs)


@pytest.fixture
def various_analysis(tmp_path):
    tmp_path = str(tmp_path)  # for Python2

    # Let's make several "reallish" test projects
    to_test = list()
    project_prefix_name = "test-project"
    data_type = "ATAC-seq"
    genome_assemblies = [("human", "hg38"), ("mouse", "mm10")]  # ("human", "hg19")
    factors = [1]  # 2, 5
    features = [100]  # 1000, 10000
    replicates = [2]  # 5

    for organism, genome_assembly in genome_assemblies:
        for n_factors in factors:
            for n_features in features:
                for n_replicates in replicates:
                    project_name = "_".join(
                        str(x) for x in [
                            project_prefix_name,
                            data_type,
                            genome_assembly,
                            n_factors,
                            n_features,
                            n_replicates])

                    an = generate_project(
                        output_dir=tmp_path,
                        project_name=project_name,
                        organism=organism,
                        genome_assembly=genome_assembly,
                        data_type=data_type,
                        n_factors=n_factors,
                        n_replicates=n_replicates,
                        n_features=n_features)
                    an.load_data()
                    to_test.append(an)
    return to_test


@pytest.fixture
def chrom_file(tmp_path):
    from ngs_toolkit.utils import download_gzip_file
    import pandas as pd

    url = (
        "https://egg2.wustl.edu/roadmap/data/byFileType/"
        + "chromhmmSegmentations/ChmmModels/coreMarks/jointModel/"
        + "final/E002_15_coreMarks_hg38lift_dense.bed.gz"
    )
    chrom_state_file = os.path.join(
        tmp_path, "E002_15_coreMarks_hg38lift_dense.bed")
    download_gzip_file(url, chrom_state_file)

    # Test
    assert os.path.exists(chrom_state_file)
    assert os.stat(chrom_state_file).st_size > 0
    b = pd.read_csv(chrom_state_file, skiprows=1, sep="\t")
    assert b.shape == (281837, 9)
    return chrom_state_file


@pytest.fixture
def get_crispr_matrix(tmp_path):
    import pandas as pd

    url = "http://liulab.dfci.harvard.edu/Mageck/melanoma.csv.zip"
    output_file = os.path.join(tmp_path, "Mageck.melanoma.example_data.zip")
    b = pd.read_csv(url, index_col=[0, 1])
    assert b.shape == (64076, 9)
    b.to_csv(output_file)
    return output_file
