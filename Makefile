.PHONY: help
help:
	@echo "help             print this message"
	@echo "check            run all tests and validations"
	@echo "check-fast       run tests and validations that finish quickly"
	@echo "check-format     verify that the code formatter has been run"
	@echo "check-slow       run tests and validations that take awhile"
	@echo "container-build  build the testing container"
	@echo "container-push   push the testing container to Docker Hub"
	@echo "format           run the code formatter"
	@echo "setup            install runtime dependencies"
	@echo "setup-dev        install development dependencies"

ifndef CONTAINER_VERSION
override CONTAINER_VERSION := latest
endif

.PHONY: check
check: check-fast check-slow check-format

.PHONY: check-fast
check-fast:
	python3 -m poetry run python -m pytest -m 'not slow' tests/ polyA/ test_inputs/AdjudicateRegions.py

.PHONY: check-format
check-format:
	python3 -m poetry run python -m black --check -t py38 -l 80 polyA/ tests/

.PHONY: check-slow
check-slow:
	python3 -m poetry run python -m pytest -m slow tests/ polyA/
	cd test_inputs && ./RunTests.sh

.PHONY: container-build
container-build:
	docker build -t traviswheelerlab/polya-build:${CONTAINER_VERSION} .

.PHONY: container-push
container-push:
	docker push traviswheelerlab/polya-build:${CONTAINER_VERSION}

.PHONY: format
format:
	python3 -m poetry run python -m black -t py38 -l 80 polyA/ tests/

.PHONY: setup
setup:
	python3 -m poetry install --no-dev

.PHONY: setup-dev
setup-dev:
	python3 -m poetry install
