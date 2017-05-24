#! /usr/bin/env python

import sys

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

version = open("VERSION").read().strip()
requirements = open("requirements.txt").read().strip().split("\n")

# setup
setup(
	name="toolkit",
	packages=["toolkit"],
	version=version,
	description="Toolkit for bioinformatics.",
	long_description=open('README.md').read(),
	classifiers=[
		"Development Status :: 3 - Alpha",
		"License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
		"Programming Language :: Python :: 2.7",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	],
	keywords="bioinformatics, sequencing, ngs, ATAC-Seq, ChIP-seq, RNA-seq",
	url="https://github.com/afrendeiro/toolkit",
	author=u"Andre Rendeiro",
	license="GPL2",
	install_requires=requirements,
	**extra
)
