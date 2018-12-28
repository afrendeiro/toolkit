.DEFAULT_GOAL := pypitest

test:
	pytest

coverage:
	coverage run -m pytest
	coverage report
	coverage html

covupload: coverage
	coverage xml
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
