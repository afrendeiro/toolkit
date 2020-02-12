#!/usr/bin/env python

import os
import sys
import subprocess
import pytest

from .conftest import CI


pytestmark = pytest.mark.skipif(
    CI, reason="Testing of recipes is not performed on CI.")


@pytest.fixture()
def pep():
    cmd = (
        "{exe} -m ngs_toolkit.recipes.generate_project"
        " --sample-input-files True").format(exe=sys.executable)
    o = subprocess.check_output(cmd.split(" "))
    o = o.decode().strip()
    assert o[-5:] == ".yaml"
    return o


def test_region_set_frip(pep):
    cmd = (
        "{exe} -m ngs_toolkit.recipes.region_set_frip {pep}"
    ).format(exe=sys.executable, pep=pep)

    p = subprocess.Popen(cmd.split(" "))
    o = p.communicate()
    assert o == (None, None)


def test_ngs_analysis(pep):
    cmd = (
        "{exe} -m ngs_toolkit.recipes.ngs_analysis {pep}"
    ).format(exe=sys.executable, pep=pep)

    p = subprocess.Popen(cmd.split(" "))
    o = p.communicate()
    assert o == (None, None)


def test_merge_signal(pep):
    from .conftest import file_exists_and_not_empty
    dir_ = os.path.dirname(os.path.dirname(pep))
    output_dir = os.path.join(dir_, "data_merged")
    cmd = (
        "{exe} -m ngs_toolkit.recipes.merge_signal "
        "--attributes A "
        "--output-dir {output_dir} "
        "{pep}"
    ).format(
        exe=sys.executable,
        output_dir=output_dir,
        pep=pep)

    p = subprocess.Popen(cmd.split(" "))
    o = p.communicate()
    assert o == (None, None)

    files = [
        "A_1.bigWig",
        "A_1.merged.bam",
        "A_1.merged.sorted.bam",
        "A_1.merged.sorted.bam.bai",
        "A_1.merge_signal.sh",
        "A_2.bigWig",
        "A_2.merged.bam",
        "A_2.merged.sorted.bam",
        "A_2.merged.sorted.bam.bai",
        "A_2.merge_signal.sh"]

    for f in files:
        file_exists_and_not_empty(os.path.join(output_dir, f))
