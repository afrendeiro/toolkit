#!/usr/bin/env python

import os

import pytest

from .conftest import CI, DEV, BUILD_DIR


def test_version_matches():
    from ngs_toolkit import __version__ as installed_version
    import pkgutil

    file_version = (
        pkgutil.get_data("ngs_toolkit", "_version.py")
        .decode()
        .strip()
        .split(" = ")[1]
        .replace('"', "")
    )
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

    replace = {"scikit-learn": "sklearn"}

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
    if CI:
        path = os.path.join(BUILD_DIR, "requirements", "requirements.txt")
    else:
        path = os.path.join("requirements", "requirements.txt")
    data = open(path).read().split("\n")

    replace = {"scikit-learn": "sklearn"}

    # handle github stuff
    requirements = [
        x.split("=")[0].replace(">", "")
        for x in data if "extra" not in x]
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
