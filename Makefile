.DEFAULT_GOAL := pypitest

build:
	python setup.py sdist bdist_wheel

pypitest: build
	twine upload -r pypitest dist/*

pypi: build
	twine upload dist/*

clean_dist:
	rm -r dist/

clean_build:
	rm -r build/

clean: clean_dist clean_build

all: build pypitest pypi clean_dist clean_build clean

.PHONY: build pypitest pypi clean_dist clean_build clean
