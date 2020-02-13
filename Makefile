.DEFAULT_GOAL := pypitest

install:
	python -m \
		pip \
		install \
		git+https://github.com/afrendeiro/toolkit.git#egg=ngs-toolkit \
		--user

test:
	python -m \
		pytest -n 3 \
		--disable-warnings \
		--show-capture=no \
		--cov=ngs_toolkit \
		--lf \
		--cov-report xml \
		ngs_toolkit/tests/test_*.py


test_cov:
	python -m \
		pytest \
		--testmon \
		--disable-warnings \
		--show-capture=no \
		ngs_toolkit/tests/test_*.py


coverage: test
	python -m codecov \
		-f coverage.xml
	python -m codacy \
		-r coverage.xml

docs:
	cd docs && $(MAKE) html
	xdg-open docs/build/html/index.html

build: test
	python setup.py sdist bdist_wheel

pypitest: build
	twine \
		upload \
		-r pypitest dist/*

pypi: build
	twine \
		upload \
		dist/*

gh:
	docker \
		build \
		-t ngs-toolkit \
		.
	docker \
		tag \
		ngs-toolkit \
		docker.pkg.github.com/afrendeiro/toolkit/ngs-toolkit:latest
	docker \
		push \
		docker.pkg.github.com/afrendeiro/toolkit/ngs-toolkit:latest

gh-release: install
	$(eval VERSION := \
		$(shell \
			python3 \
			-c 'from ngs_toolkit import __version__ as v; print(v)'))
	docker \
		build \
		-t ngs-toolkit:$(VERSION) \
		.
	docker \
		tag \
		ngs-toolkit \
		docker.pkg.github.com/afrendeiro/toolkit/ngs-toolkit:$(VERSION)
	docker \
		push \
		docker.pkg.github.com/afrendeiro/toolkit/ngs-toolkit:$(VERSION)

clean_pyc:
	find . -name \*.pyc -delete

clean_mypy:
	rm -rf .mypy_cache/

clean_test:
	rm -rf .pytest_cache/
	rm -rf /tmp/pytest*
	find . -name "__pycache__" -exec rm -rf {} \;
	rm -rf .coverage*
	rm -rf .tox/

clean_cov: clean_test
	rm -fr coverage.xml htmlcov

clean_docs:
	rm -fr docs/build/

clean_dist:
	rm -fr dist/

clean_build:
	rm -fr build/
	rm -rf ngs_toolkit/_version.py

clean_eggs:
	rm -fr ngs_toolkit.egg-info
	rm -fr .eggs

clean: \
	clean_pyc \
	clean_mypy \
	clean_test \
	clean_cov \
	clean_docs \
	clean_dist \
	clean_build \
	clean_eggs

all: \
	test \
	coverage \
	docs \
	build \
	pypitest \
	pypi \
	clean

.PHONY: \
	test \
	coverage \
	docs \
	build \
	pypitest \
	pypi \
	clean_pyc \
	clean_mypy \
	clean_test \
	clean_cov \
	clean_docs \
	clean_dist \
	clean_build \
	clean_eggs \
	clean
