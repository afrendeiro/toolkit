#!/usr/bin/env python

import os
import sys
import subprocess
import pytest
import pandas as pd

from ngs_toolkit import Analysis
from .conftest import CI, DEV, file_exists_and_not_empty


def test_region_set_frip(pep):
    import pkgutil

    # For this test, we need an analysis object with sample attributes pointing
    # to their input files (just like the atac_analysis_with_input_files parent),
    # but since it has to be in inside the recipe, so we will temporarily set it
    # at the home directory level, for this test only
    config = os.path.join(os.path.expanduser("~"), ".ngs_toolkit.config.yaml")
    yaml = pkgutil.get_data("ngs_toolkit", "config/example.yaml").decode().strip()
    with open(config, "w") as handle:
        handle.write(yaml)

    cmd = (
        "{exe} -m ngs_toolkit.recipes.region_set_frip {pep} " "--computing-configuration default"
    ).format(exe=sys.executable, pep=pep)

    p = subprocess.Popen(cmd.split(" "))
    o = p.communicate()
    assert o == (None, None)

    an = Analysis(from_pep=pep)
    for sample in an.samples:
        files = [
            "region_set_frip.all_reads.txt",
            "region_set_frip.inside_reads.txt",
            sample.name + ".region_set_frip.log",
            sample.name + ".region_set_frip.sh",
            "stats.tsv",
        ]
        for f in files:
            assert file_exists_and_not_empty(os.path.join(sample.sample_root, f))

    os.remove(config)


def test_deseq2(tmp_path, atac_analysis_with_input_files):
    an = atac_analysis_with_input_files
    an.differential_analysis(distributed=True, dry_run=True)

    p = "differential_analysis_ATAC-seq"
    comp_dir = os.path.join(an.results_dir, p, "Factor_A_2vs1")
    output_prefix = os.path.join(comp_dir, p)

    cmd = ("{} -m ngs_toolkit.recipes.deseq2 --output-prefix {} {}").format(
        sys.executable, output_prefix, comp_dir
    )

    proc = subprocess.Popen(cmd.split(" "))
    o = proc.communicate()
    assert o == (None, None)

    files = [
        # "deseq_job.Factor_A_2vs1.log",
        p + ".deseq_result.Factor_A_2vs1.csv",
        p + ".deseq_result.all_comparisons.csv",
    ]
    for f in files:
        assert file_exists_and_not_empty(os.path.join(comp_dir, f))


def test_coverage(tmp_path, atac_analysis_with_input_files):

    region_set = atac_analysis_with_input_files.sites.fn
    sample = atac_analysis_with_input_files.samples[0]
    output = os.path.join(tmp_path, "output.bed")

    cmd = ("{} -m ngs_toolkit.recipes.coverage {} {} {}").format(
        sys.executable, region_set, sample.aligned_filtered_bam, output
    )

    p = subprocess.Popen(cmd.split(" "))
    o = p.communicate()
    assert o == (None, None)

    assert file_exists_and_not_empty(output)

    assert pd.read_csv(output, sep="\t", header=None).shape[1] == 4


def test_enrichr_good(tmp_path):
    genes = ["PAX5", "SOX2"]
    input_file = os.path.join(tmp_path, "genes.txt")
    output_file = os.path.join(tmp_path, "enrichr.csv")
    with open(input_file, "w") as handle:
        for g in genes:
            handle.write(g + "\n")
    cmd = ("{} -m ngs_toolkit.recipes.enrichr {} {}").format(
        sys.executable, input_file, output_file
    )

    p = subprocess.Popen(cmd.split(" "))
    o = p.communicate()
    assert o == (None, None)

    assert file_exists_and_not_empty(output_file)

    assert pd.read_csv(output_file).shape[1] == 10


def test_enrichr_bad(tmp_path):
    # No genes enriched, should return empty dataframe
    genes = ["!!~~IMPOSSIBLEGENE~~!!"]
    input_file = os.path.join(tmp_path, "impossible_genes.txt")
    output_file = os.path.join(tmp_path, "empty_enrichr.csv")
    with open(input_file, "w") as handle:
        for g in genes:
            handle.write(g + "\n")
    cmd = ("{} -m ngs_toolkit.recipes.enrichr {} {}").format(
        sys.executable, input_file, output_file
    )

    p = subprocess.Popen(cmd.split(" "))
    o = p.communicate()
    assert o == (None, None)

    assert file_exists_and_not_empty(output_file)

    with pytest.raises(pd.errors.EmptyDataError):
        pd.read_csv(output_file)


@pytest.mark.skipif(CI or DEV, reason="Test too long to be performed on CI.")
def test_ngs_analysis(pep):
    cmd = ("{exe} -m ngs_toolkit.recipes.ngs_analysis {pep}").format(exe=sys.executable, pep=pep)

    p = subprocess.Popen(cmd.split(" "))
    o = p.communicate()
    assert o == (None, None)


@pytest.mark.skipif(CI or DEV, reason="Test too long to be performed on CI.")
def test_merge_signal(pep):
    import pkgutil
    from .conftest import file_exists_and_not_empty

    dir_ = os.path.dirname(os.path.dirname(pep))
    output_dir = os.path.join(dir_, "data_merged")
    cmd = (
        "{exe} -m ngs_toolkit.recipes.merge_signal "
        "-d "
        "--attributes A "
        "--output-dir {output_dir} "
        "{pep}"
    ).format(exe=sys.executable, output_dir=output_dir, pep=pep)

    # this requires a config with sample input files
    file_config = os.path.join(os.path.expanduser("~"), ".ngs_toolkit.config.yaml")
    content = pkgutil.get_data("ngs_toolkit", "config/default.yaml").decode().strip()
    with open(file_config, "w") as handle:
        handle.write(content)

    p = subprocess.Popen(cmd.split(" "))
    o = p.communicate()
    assert o == (None, None)

    files = [
        # "A_1.bigWig",
        # "A_1.merged.bam",
        # "A_1.merged.sorted.bam",
        # "A_1.merged.sorted.bam.bai",
        "A_1.merge_signal.sh",
        # "A_2.bigWig",
        # "A_2.merged.bam",
        # "A_2.merged.sorted.bam",
        # "A_2.merged.sorted.bam.bai",
        "A_2.merge_signal.sh",
    ]

    for f in files:
        assert file_exists_and_not_empty(os.path.join(output_dir, f))

    os.remove(file_config)
