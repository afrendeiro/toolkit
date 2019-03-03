#!/usr/bin/env python


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


def test_all_requirements_are_importable():
    import requests
    import importlib

    package_name = "ngs-toolkit"
    url = "https://pypi.python.org/pypi/" + str(package_name) + "/json"
    data = requests.get(url).json()

    replace = {"scikit-learn": "sklearn", "piper": "pypiper"}

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
