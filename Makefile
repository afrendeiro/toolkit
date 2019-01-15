.DEFAULT_GOAL := pypitest

test3:
	python3 -m pytest -n 3 --disable-warnings --show-capture=no --cov=ngs_toolkit --cov-report xml test/test_*.py --lf

test2:
	python2 -m pytest -n 3 --disable-warnings --show-capture=no --cov=ngs_toolkit --cov-report xml test/test_*.py --lf

test: test3

coverage: test
	codecov -f coverage.xml
	python-codacy-coverage -r coverage.xml

build: test
	python setup.py sdist bdist_wheel

pypitest: build
	twine upload -r pypitest dist/*

pypi: build
	twine upload dist/*

clean_cov:
	rm -r coverage.xml htmlcov

clean_dist:
	rm -r dist/

clean_build:
	rm -r build/

clean: clean_cov clean_dist clean_build

all: test coverage build pypitest pypi clean_dist clean_build clean

.PHONY: test coverage build pypitest pypi clean_dist clean_build clean
