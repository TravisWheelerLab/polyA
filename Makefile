.PHONY: help
help:
	@echo "help             print this message"
	@echo "check            run all tests and validations"
	@echo "check-fast       run tests and validations that finish quickly"
	@echo "check-format     verify that the code formatter has been run"
	@echo "check-slow       run tests and validations that take awhile"
	@echo "container        build and push a new container image"
	@echo "container-build  build the testing container"
	@echo "container-push   push the testing container to Docker Hub"
	@echo "format           run the code formatter"
	@echo "setup            install runtime dependencies"
	@echo "setup-dev        install development dependencies"

ifndef CONTAINER_VERSION
override CONTAINER_VERSION := latest
endif

RUN_CMD := pipenv run
PYTHON_CMD := ${RUN_CMD} python

FMT_CMD := ${PYTHON_CMD} -m black
FMT_TARGETS := polyA/ tests/
FMT_OPTS := -t py38 -l 80

TEST_CMD := ${PYTHON_CMD} -m pytest
TEST_TARGETS := tests/ polyA/

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
check: check-fast check-slow check-format check-types

.PHONY: check-benchmark
check-benchmark:
	cd tests && PYTHONPATH=../ ${RUN_CMD} ./RunBenchmarks.sh

.PHONY: check-fast
check-fast:
	PYTHONPATH=./ ${TEST_CMD} ${TEST_TARGETS}

.PHONY: check-format
check-format:
	${FMT_CMD} --check ${FMT_OPTS} ${FMT_TARGETS}

.PHONY: check-slow
check-slow:
	cd tests && PYTHONPATH=../ ${RUN_CMD} ./RunTests.sh
	cd tests && PYTHONPATH=../ ${RUN_CMD} ./RunUltraTests.sh

.PHONY: check-types
check-types:
	${PYTHON_CMD} -m mypy -m polyA

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

.PHONY: format
format:
	${FMT_CMD} ${FMT_OPTS} ${FMT_TARGETS}

.PHONY: setup
setup:
	pipenv install

.PHONY: setup-dev
setup-dev:
	pipenv install --dev
