#!/usr/bin/env python

import os

import pytest

from .conftest import CI, DEV, BUILD_DIR


@pytest.mark.skipif(
    not CI,
    reason="Development mode, not testing Pypi requirements")
def test_version_matches():
    from ngs_toolkit import __version__ as installed_version
    from pkg_resources import get_distribution

    file_version = get_distribution('ngs_toolkit').version

    assert installed_version == file_version


@pytest.mark.skipif(
    DEV,
    reason="Development mode, not testing Pypi requirements")
def test_pypi_requirements_are_importable():
    import requests
    import importlib

    package_name = "ngs-toolkit"
    url = "https://pypi.python.org/pypi/" + str(package_name) + "/json"
    data = requests.get(url).json()

    # handle packages where the package name is different than the Pypi name
    replace = {
        "setuptools-scm": "setuptools_scm",
        "scikit-learn": "sklearn"}

    # not extra requirements
    requirements = [
        x.split(" ")[0]
        for x in data["info"]["requires_dist"]
        if "extra" not in x
    ]
    for req in requirements:
        if req in replace:
            requirements.pop(requirements.index(req))
            requirements.append(replace[req])

    for req in requirements:
        try:
            importlib.import_module(req)
        except ImportError:
            assert False


def test_all_requirements_are_importable():
    import importlib

    # test only basic requirements (not extras)
    path = None
    if CI:
        reqs = os.path.join(BUILD_DIR, "requirements", "requirements.txt")
        if os.path.exists(reqs):
            path = reqs
    if path is None:
        path = os.path.join("requirements", "requirements.txt")

    if not os.path.exists(path):
        pytest.skip("Could not locate requirements.txt")

    data = open(path).read().split("\n")

    replace = {"scikit-learn": "sklearn"}

    # handle github stuff
    requirements = list()
    for x in data:
        for sep in ['>=', '<=', '=', ">", "<"]:
            x = x.split(sep)
            if len(x) == 2:
                x = x[0].replace("=", "").replace(">", "").replace("<", "")
            else:
                x = x[0]
        if "extra" not in x:
            requirements.append(x)

    # remove commnets
    requirements = [x[:x.index("  #")] if "#" in x else x
                    for x in requirements]
    # remove empty lines
    requirements = [x for x in requirements if x != ""]
    for req in requirements:
        if req in replace:
            requirements.pop(requirements.index(req))
            requirements.append(replace[req])

    for req in requirements:
        try:
            importlib.import_module(req)
        except ImportError:
            assert False, "Required '%s' module could not be found!" % req
