.DEFAULT_GOAL := pypitest

install:
	pip install .[testing] --user

test:
	python -m pytest -n 3 --disable-warnings --show-capture=no --cov=ngs_toolkit --cov-report xml ngs_toolkit/tests/test_*.py --lf

coverage: test
	codecov -f coverage.xml
	python-codacy-coverage -r coverage.xml

docs:
	cd docs && $(MAKE) html
	xdg-open docs/build/html/index.html

build: test
	python setup.py sdist bdist_wheel

pypitest: build
	twine upload -r pypitest dist/*

pypi: build
	twine upload dist/*

clean_test:
	rm -rf .pytest_cache/
	find . -name "__pycache__" -exec rm -r {} \;
	rm .coverage*

clean_cov: clean_test
	rm -fr coverage.xml htmlcov

clean_docs:
	rm -fr docs/build/

clean_dist:
	rm -fr dist/

clean_build:
	rm -fr build/
	rm -fr ngs_toolkit.egg-info
	find . -name \*.pyc -delete

clean: clean_test clean_cov clean_docs clean_dist clean_build

all: test coverage docs build pypitest pypi clean_test clean_cov clean_docs clean_dist clean_build clean

.PHONY: test coverage docs build pypitest pypi clean_test clean_cov clean_docs clean_dist clean_build clean
