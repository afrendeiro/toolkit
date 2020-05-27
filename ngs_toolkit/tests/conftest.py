#!/usr/bin/env python


import os
from functools import partialmethod

import pytest

from ngs_toolkit import MEMORY, _CONFIG, _LOGGER, Analysis, ATACSeqAnalysis
from ngs_toolkit.demo.data_generator import generate_project


# Environment-specific
CI: bool = ("TRAVIS" in os.environ) or ("GITHUB_WORKFLOW" in os.environ)
CI_NAME = None
BUILD_DIR: str = os.path.abspath(os.path.curdir)
DEV: bool = False
RPY2: bool
DESEQ2: bool
PREPROCESSCORE: bool
R: bool
COMBAT: bool


if CI:
    if "TRAVIS" in os.environ:
        CI_NAME = "TRAVIS"
        BUILD_DIR = os.environ["TRAVIS_BUILD_DIR"]
    elif "GITHUB_WORKFLOW" in os.environ:
        CI_NAME = "GITHUB"
        BUILD_DIR = os.path.join("home", "runner", "work", "toolkit", "toolkit")

try:
    DEV = os.environ["TRAVIS_BRANCH"] == "dev"
except KeyError:
    pass
try:
    DEV = os.environ["GITHUB_REF"] == "dev"
except KeyError:
    import subprocess

    try:
        o = subprocess.check_output("git status".split(" "))
        DEV = "dev" in o.decode().split("\n")[0]
    except subprocess.CalledProcessError:
        msg = "Could not detect whether on a development branch."
        _LOGGER.warning(msg)


# Test-specifc options
# # Note:
# # The DESeq2 1.24.0 version in Debian archives
# # differs from the DESeq2 1.24.0 version in bioconductor version 3.9
# # If estimateDispersions with default fitType="parametric" fails,
# # (as often happens with the quickly generated synthetic data from tests),
# # it tries to use local fit using the locfit package, but in Debian
# # version this is not a valid choice of fit, causing failure.
# # Due to this, and since I'm using Debian packages for faster testing
# # I'm manually setting fitType="mean" for testing only.
Analysis.differential_analysis = partialmethod(
    Analysis.differential_analysis, deseq_kwargs={"fitType": "mean"}
)


# This a part of the "example" config that is required for some analysis
# that require existing input files (BAM, peaks, etc)
NEW_CONFIG = {
    "sample_input_files": {
        "ATAC-seq": {
            "aligned_filtered_bam": "{data_dir}/{sample_name}/mapped/{sample_name}.trimmed.bowtie2.filtered.bam",
            "peaks": "{data_dir}/{sample_name}/peaks/{sample_name}_peaks.narrowPeak",
            "summits": "{data_dir}/{sample_name}/peaks/{sample_name}_summits.bed",
        },
        "ChIP-seq": {
            "aligned_filtered_bam": "{data_dir}/{sample_name}/mapped/{sample_name}.trimmed.bowtie2.filtered.bam"
        },
        "CNV": {
            "log2_read_counts": {
                "10kb": "{data_dir}/{sample_name}/{sample_name}_10kb/CNAprofiles/log2_read_counts.igv",
                "100kb": "{data_dir}/{sample_name}/{sample_name}_100kb/CNAprofiles/log2_read_counts.igv",
                "1000kb": "{data_dir}/{sample_name}/{sample_name}_1000kb/CNAprofiles/log2_read_counts.igv",
            }
        },
        "RNA-seq": {
            "aligned_filtered_bam": "{data_dir}/{sample_name}/mapped/{sample_name}.trimmed.bowtie2.filtered.bam",
            "bitseq_counts": "{data_dir}/{sample_name}/bowtie1_{genome}/bitSeq/{sample_name}.counts",
        },
    }
}


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


def has_module(module):
    import importlib

    try:
        importlib.import_module(module)
        return True
    except ModuleNotFoundError:
        return False


def has_R_library(library):
    if not has_module("rpy2"):
        return False
    from rpy2.robjects.packages import importr
    from rpy2.robjects.packages import PackageNotInstalledError

    try:
        importr(library)
        return True
    except PackageNotInstalledError:
        return False


RPY2 = has_module("rpy2")
COMBAT = has_module("combat")
STAP = has_module("STAP")
DNACOPY = has_module("DNAcopy")
PREPROCESSCORE = has_R_library("preprocessCore")
DESEQ2 = has_R_library("DESeq2")
R = RPY2 and PREPROCESSCORE and DESEQ2
R_REASON = "R, rpy2 or R libraries not available."


def file_exists(file):
    from ngs_toolkit.utils import get_this_file_or_timestamped

    return os.path.exists(get_this_file_or_timestamped(file))


def file_not_empty(file):
    from ngs_toolkit.utils import get_this_file_or_timestamped

    return os.stat(get_this_file_or_timestamped(file)).st_size > 0


def file_exists_and_not_empty(file):
    from ngs_toolkit.utils import get_this_file_or_timestamped

    f = get_this_file_or_timestamped(file)

    return os.path.exists(f) and (os.stat(f).st_size > 0)


@pytest.fixture
def empty_analysis():
    return ATACSeqAnalysis()


@pytest.fixture
def null_analysis():
    an = ATACSeqAnalysis()
    an.organism = None
    an.genome = None
    an.sites = None
    return an


@pytest.fixture
def full_analysis():
    an = ATACSeqAnalysis()
    an.organism = "human"
    an.genome = "hg38"
    an.sites = "hg38"
    an.samples = []
    return an


@pytest.fixture
def analysis():
    return Analysis(name="test-project")


@pytest.fixture
def atac_analysis(tmp_path):
    tmp_path = str(tmp_path)
    kwargs = {
        "data_type": "ATAC-seq",
        "organism": "human",
        "genome_assembly": "hg38",
        "n_factors": 1,
        "n_features": 250,
        "n_replicates": 2,
    }
    kwargs.update(
        {
            "project_name": "test-project_" + "_".join(str(x) for x in kwargs.values()),
            "output_dir": tmp_path,
        }
    )

    return generate_project(**kwargs)


@pytest.fixture
def atac_analysis_with_input_files(tmp_path):
    from ngs_toolkit.demo.data_generator import generate_sample_input_files

    _CONFIG.update(NEW_CONFIG)

    tmp_path = str(tmp_path)
    kwargs = {
        "data_type": "ATAC-seq",
        "organism": "human",
        "genome_assembly": "hg38",
        "n_factors": 1,
        "n_features": 10,
        "n_replicates": 2,
    }
    kwargs.update(
        {
            "project_name": "test-project_" + "_".join(str(x) for x in kwargs.values()),
            "output_dir": tmp_path,
        }
    )
    an = generate_project(**kwargs)
    an.set_samples_input_files()
    generate_sample_input_files(an, an.matrix_raw)

    return an


@pytest.fixture
def atac_analysis_with_unmapped_input_files(atac_analysis_with_input_files):
    import pandas as pd
    import yaml

    with atac_analysis_with_input_files as an:
        # We will update the annotation to add a 'data_source' column
        csv = an.pep.replace("project_config.yaml", "annotation.csv")
        df = pd.read_csv(csv)
        df["data_source"] = "mapped"
        df.to_csv(csv, index=False)
        # We will update the config to add a line pointing to the aligned bams
        conf = yaml.safe_load(open(an.pep))
        mapped_path = os.path.join(
            an.root_dir, "data/{sample_name}/mapped/{sample_name}.trimmed.bowtie2.filtered.bam"
        )
        conf["sample_modifiers"]["derive"]["attributes"].append("mapped")
        conf["sample_modifiers"]["derive"]["sources"]["mapped"] = mapped_path

        yaml.safe_dump(conf, open(an.pep, "w"))
        a = ATACSeqAnalysis(from_pep=an.pep)
        a.load_data()

        return a


@pytest.fixture
def subproject_config(atac_analysis):
    import yaml

    annot = os.path.join(atac_analysis.root_dir, "metadata", "annotation.csv")
    subannot = os.path.join(atac_analysis.root_dir, "metadata", "sample_subannotation.csv")

    yaml_file = os.path.join(atac_analysis.root_dir, "metadata", "project_config.yaml")
    conf = yaml.safe_load(open(yaml_file, "r"))
    conf["project_modifiers"] = {
        "amend": {"test_subproject": {"sample_table": annot, "subsample_table": subannot}}
    }
    del conf["sample_table"]
    del conf["subsample_table"]

    yaml.safe_dump(conf, open(yaml_file, "w"))

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
    tmp_path = str(tmp_path)
    kwargs = {
        "data_type": "ATAC-seq",
        "organism": "human",
        "genome_assembly": "hg38",
        "n_factors": 4,
        "n_features": 100,
        "n_replicates": 4,
    }
    kwargs.update(
        {
            "project_name": "test-project_" + "_".join(str(x) for x in kwargs.values()),
            "output_dir": tmp_path,
        }
    )

    an = generate_project(**kwargs)
    an.load_data()
    an.normalize(method="rpm")
    an.annotate_samples()

    return an


@pytest.fixture
def rnaseq_analysis(tmp_path):
    tmp_path = str(tmp_path)
    kwargs = {
        "data_type": "RNA-seq",
        "organism": "human",
        "genome_assembly": "hg38",
        "n_factors": 1,
        "n_features": 100,
        "n_replicates": 3,
    }
    kwargs.update(
        {
            "project_name": "test-project_" + "_".join(str(x) for x in kwargs.values()),
            "output_dir": tmp_path,
        }
    )

    return generate_project(**kwargs)


@pytest.fixture
def chipseq_analysis(tmp_path):
    tmp_path = str(tmp_path)

    _CONFIG.update(NEW_CONFIG)

    kwargs = {
        "data_type": "ChIP-seq",
        "organism": "human",
        "genome_assembly": "hg38",
        "n_factors": 2,
        "n_features": 100,
        "n_replicates": 2,
    }
    kwargs.update(
        {
            "project_name": "test-project_" + "_".join(str(x) for x in kwargs.values()),
            "output_dir": tmp_path,
        }
    )

    an = generate_project(**kwargs, sample_input_files=True)
    # an = generate_project(**kwargs, sample_input_files=True)
    an.comparison_table["comparison_type"] = "peaks"
    an.set_comparisons()

    return an


@pytest.fixture
def chipseq_analysis_with_peaks(chipseq_analysis):
    import numpy as np
    from ngs_toolkit.utils import bed_to_index

    df = chipseq_analysis.sites.to_dataframe()

    # for homer peaks add p-value column
    df["name"] = bed_to_index(df)  # dummy column
    df["score"] = np.random.random(df.shape[0])  # p-value
    df["strand"] = "."  # dummy column

    for name, comp in chipseq_analysis.comparisons.items():

        os.makedirs(comp["output_dir"])

        for peak_type in ["original", "filtered"]:
            for peak_caller, file in comp["peak_calls"][peak_type].items():
                # select 30% of all peaks to be present in each sample
                df2 = df.sample(frac=0.3)

                # for homer the column order is different
                if "homer" in peak_caller:
                    df2 = df2[["name", "chrom", "start", "end", "score", "strand"]]

                with open(file, "w") as handle:
                    for _, entry in df2.iterrows():
                        handle.write("\t".join(entry.astype(str)) + "\n")

    return chipseq_analysis


@pytest.fixture
def cnv_analysis(tmp_path):
    tmp_path = str(tmp_path)

    _CONFIG.update(NEW_CONFIG)

    kwargs = {
        "data_type": "CNV",
        "organism": "human",
        "genome_assembly": "hg38",
        "n_factors": 2,
        "n_features": 100,
        "n_replicates": 2,
    }
    kwargs.update(
        {
            "project_name": "test-project_" + "_".join(str(x) for x in kwargs.values()),
            "output_dir": tmp_path,
        }
    )

    an = generate_project(**kwargs, sample_input_files=False)
    return an


@pytest.fixture
def cnv_analysis_with_inputs(tmp_path):
    tmp_path = str(tmp_path)

    _CONFIG.update(NEW_CONFIG)

    kwargs = {
        "data_type": "CNV",
        "organism": "human",
        "genome_assembly": "hg38",
        "n_factors": 2,
        "n_features": 100,
        "n_replicates": 2,
    }
    kwargs.update(
        {
            "project_name": "test-project_" + "_".join(str(x) for x in kwargs.values()),
            "output_dir": tmp_path,
        }
    )

    an = generate_project(**kwargs, sample_input_files=True)
    return an


@pytest.fixture()
def pep(atac_analysis_with_input_files):
    return atac_analysis_with_input_files.pep


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
                        str(x)
                        for x in [
                            project_prefix_name,
                            data_type,
                            genome_assembly,
                            n_factors,
                            n_features,
                            n_replicates,
                        ]
                    )

                    an = generate_project(
                        output_dir=tmp_path,
                        project_name=project_name,
                        organism=organism,
                        genome_assembly=genome_assembly,
                        data_type=data_type,
                        n_factors=n_factors,
                        n_replicates=n_replicates,
                        n_features=n_features,
                    )
                    an.load_data()
                    to_test.append(an)
    return to_test


@pytest.fixture
@MEMORY.cache
def chrom_file(tmp_path):
    from ngs_toolkit.utils import download_gzip_file
    import pandas as pd

    url = (
        "https://egg2.wustl.edu/roadmap/data/byFileType/"
        + "chromhmmSegmentations/ChmmModels/coreMarks/jointModel/"
        + "final/E002_15_coreMarks_hg38lift_dense.bed.gz"
    )
    chrom_state_file = os.path.join(tmp_path, "E002_15_coreMarks_hg38lift_dense.bed")
    download_gzip_file(url, chrom_state_file)

    # Test
    assert os.path.exists(chrom_state_file)
    assert os.stat(chrom_state_file).st_size > 0
    b = pd.read_csv(chrom_state_file, skiprows=1, sep="\t")
    assert b.shape == (281837, 9)
    return chrom_state_file


@pytest.fixture
@MEMORY.cache
def get_crispr_matrix(tmp_path):
    import pandas as pd

    url = (
        "https://web.archive.org/web/20190702030135/"
        "http://liulab.dfci.harvard.edu/Mageck/melanoma.csv.zip"
    )
    b = pd.read_csv(url, index_col=[0, 1])
    assert b.shape == (64076, 9)
    output_file = os.path.join(tmp_path, "Mageck.melanoma.example_data.zip")
    b.to_csv(output_file)
    return output_file
