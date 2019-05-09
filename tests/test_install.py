#!/usr/bin/env python

import os

import pytest


travis = "TRAVIS" in os.environ
if travis:
    dev = os.environ['TRAVIS_BRANCH'] == 'dev'
else:
    import subprocess
    o = subprocess.check_output("git status".split(" "))
    dev = "dev" in o.decode().split("\n")[0]


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


@pytest.mark.skipif(dev, reason="Development mode, not testing Pypi requirements")
def test_pypi_requirements_are_importable():
    import requests
    import importlib

    package_name = "ngs-toolkit"
    url = "https://pypi.python.org/pypi/" + str(package_name) + "/json"
    data = requests.get(url).json()

    replace = {"scikit-learn": "sklearn"}

    # not extra requirements
    requirements = [
        x.split(" ")[0] for x in data["info"]["requires_dist"] if "extra" not in x
    ]
    # versions = [x.split(" ")[1] if len(x.split(" ")) > 1 else "" for x in data['info']['requires_dist']]

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

    # not extra requirements
    data = open("requirements/requirements.txt").read().split()

    replace = {"scikit-learn": "sklearn"}

    requirements = [
        x.split("=")[0].replace(">", "") for x in data if "extra" not in x
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
