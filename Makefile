.PHONY: help
help:
	@echo "help             print this message"
	@echo "build-package    build a PyPI package"
	@echo "publish-package  publish the package to PyPI"
	@echo "check            run all tests and validations"
	@echo "check-fast       run tests and validations that finish quickly"
	@echo "check-format     verify that the code formatter has been run"
	@echo "check-types      run mypy to check type annotations"
	@echo "check-slow       run tests and validations that take awhile"
	@echo "container        build and push a new container image"
	@echo "container-build  build the testing container"
	@echo "container-push   push the testing container to Docker Hub"
	@echo "docs             build HTML version of documentation"
	@echo "docs-serve       serve the HTML documentation on port 8000"
	@echo "format           run the code formatter"
	@echo "setup            install runtime dependencies"
	@echo "setup-dev        install development dependencies"

ifndef CONTAINER_VERSION
override CONTAINER_VERSION := latest
endif

RUN_CMD := pipenv run
PYTHON_CMD := ${RUN_CMD} python

DOCS_CMD := ${RUN_CMD} sphinx-apidoc
DOCS_OPTS := -f -o docs/source polyA

FMT_CMD := ${PYTHON_CMD} -m black
FMT_TARGETS := polyA/ tests/
FMT_OPTS := -t py38 -l 80

TEST_CMD := ${PYTHON_CMD} -m pytest
TEST_TARGETS := tests/ polyA/

.PHONY: containerized
containerized:
	docker run --mount src="${PWD}",target=/code,type=bind traviswheelerlab/polya-build ${TASK}

.PHONY: build-package
build-package:
	pipenv run flit build

.PHONY: publish-package
publish-package:
	pipenv run flit publish

.PHONY: build-conda-package
build-conda-package:
	conda-build .
	mkdir -p conda-dist/
	conda convert --platform all "$(shell conda-build --output .)" -o conda-dist/

.PHONY: publish-conda-package
publish-conda-package:
	@echo "NOT IMPLEMENTED"

.PHONY: check
check: check-format check-lints check-types check-fast check-slow

.PHONY: check-fast
check-fast:
	PYTHONPATH=./ ${TEST_CMD} ${TEST_TARGETS}

.PHONY: check-format
check-format:
	${FMT_CMD} --check ${FMT_OPTS} ${FMT_TARGETS}

.PHONY: check-lints
check-lints:
	${RUN_CMD} pylint -E polyA/

.PHONY: check-types
check-types:
	${RUN_CMD} mypy polyA/

.PHONY: check-slow
check-slow:
	cd tests && PYTHONPATH=../ ${RUN_CMD} ./RunTests.sh
	cd tests && PYTHONPATH=../ ${RUN_CMD} ./RunUltraTests.sh

.PHONY: clean
clean:
	rm -rf dist/

.PHONY: container
container: container-build container-push

.PHONY: container-build
container-build:
	docker build -t traviswheelerlab/polya-build:${CONTAINER_VERSION} .

.PHONY: container-push
container-push:
	docker push traviswheelerlab/polya-build:${CONTAINER_VERSION}

.PHONY: container-build-conda-package
container-build-conda-package:
	docker run --mount src="${PWD}",target=/code,type=bind traviswheelerlab/polya-conda build-conda-package

.PHONY: container-conda
container-conda: container-conda-build container-conda-push

.PHONY: container-conda-build
container-conda-build:
	docker build -t traviswheelerlab/polya-conda:${CONTAINER_VERSION} -f Dockerfile_conda .

.PHONY: container-conda-push
container-conda-push:
	docker push traviswheelerlab/polya-conda:${CONTAINER_VERSION}

.PHONY: container-esl_scorematrix
container-esl_scorematrix: container-esl_scorematrix-build container-esl_scorematrix-push

.PHONY: container-esl_scorematrix-build
container-esl_scorematrix-build:
	docker build -t traviswheelerlab/polya-esl_scorematrix:${CONTAINER_VERSION} -f Dockerfile_esl_scorematrix .

.PHONY: container-esl_scorematrix-push
container-esl_scorematrix-push:
	docker push traviswheelerlab/polya-esl_scorematrix:${CONTAINER_VERSION}

.PHONY: docs
docs:
	${DOCS_CMD} ${DOCS_OPTS}
	${RUN_CMD} pip freeze > requirements.txt
	cd docs && make html

.PHONY: docs-serve
docs-serve:
	${PYTHON_CMD} -m http.server --directory docs/build/html

.PHONY: format
format:
	${FMT_CMD} ${FMT_OPTS} ${FMT_TARGETS}

.PHONY: setup
setup:
	pipenv install

.PHONY: setup-dev
setup-dev:
	pipenv install --dev
