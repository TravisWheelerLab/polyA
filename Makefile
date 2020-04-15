.PHONY: help
help:
	@echo "help           print this message"
	@echo "check          run all tests and validations"
	@echo "check-fast     run tests and validations that finish quickly"
	@echo "check-format   verify that the code formatter has been run"
	@echo "check-slow     run tests and validations that take awhile"
	@echo "container      build and push the CI container"
	@echo "format         run the code formatter"
	@echo "setup          install runtime dependencies"
	@echo "setup-dev      install development dependencies"

.PHONY: check
check: check-fast check-slow check-format

.PHONY: check-fast
check-fast:
	poetry run python -m pytest -m 'not slow' tests/ polyA/

.PHONY: check-format
check-format:
	poetry run python -m black --check .

.PHONY: check-slow
check-slow:
	poetry run python -m pytest -m slow tests/ polyA/
	cd test_inputs && ./RunTests.sh

.PHONY: container
container:
	docker build -t traviswheelerlab/polya-build:latest .
	docker push traviswheelerlab/polya-build:latest

.PHONY: format
format:
	poetry run python -m black .

.PHONY: setup
setup:
	poetry install --no-dev

.PHONY: setup-dev
setup-dev:
	poetry install
