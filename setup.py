#! /usr/bin/env python

import sys
import os
import glob


# take care of extra required modules depending on Python version
extra = {}

try:
    from setuptools import setup
    if sys.version_info < (2, 7):
        extra['install_requires'] = ['argparse']
    if sys.version_info >= (3,):
        extra['use_2to3'] = True
except ImportError:
    from distutils.core import setup
    if sys.version_info < (2, 7):
        extra['dependencies'] = ['argparse']

with open(os.path.join("ngs_toolkit", "_version.py"), 'r') as handle:
    version = handle.readline().split()[-1].strip("\"'\n")

requirements = open("requirements.txt").read().strip().split("\n")
requirements = [r for r in requirements if not r.startswith("#")]
requirements = [r for r in requirements if "#egg=" not in r]


# Recipes
recipes = glob.glob(os.path.join(os.path.curdir, "ngs_toolkit", "recipes", "*.py"))
recipes = list(map(
    lambda x: x.replace(".py", "").replace("./", "").replace("/", "."),
    recipes))
recipes = [" = ".join([i, j + ":main"]) for i, j in zip(map(lambda x: x.split('.')[-1], recipes), recipes)]
recipes = [r for r in recipes if not r.startswith("__init__")]


# Handle the pypi README formatting
try:
    import pypandoc
    long_description = pypandoc.convert_file('README.md', 'rst')
except(IOError, ImportError):
    long_description = open('README.md').read()


# setup
setup(
    name="ngs_toolkit",
    packages=["ngs_toolkit", "ngs_toolkit.recipes"],
    version=version,
    entry_points={
        "console_scripts": [
            'projectmanager = ngs_toolkit.project_manager:main',
            'trackmanager = ngs_toolkit.track_manager:main'
        ] + recipes,
    },
    description="A toolkit for NGS analysis with Python.",
    long_description=long_description,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    keywords="bioinformatics, sequencing, ngs, ngs analysis, ATAC-Seq, ChIP-seq, RNA-seq, project management",
    url="https://github.com/afrendeiro/toolkit",
    author=u"Andre Rendeiro",
    author_email='afrendeiro@gmail.com',
    license="GPL3",
    install_requires=requirements,
    extras_require={
        'deseq2':  ["rpy2==2.8.6"]},
    package_data={
        'ngs_toolkit': ['config/*.yaml']
    },
    **extra
)
