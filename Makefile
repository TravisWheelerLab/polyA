.PHONY: help
help:
	@echo "help           print this message"
	@echo "check          run all tests and validations"
	@echo "check-fast     run tests and validations that finish quickly"
	@echo "check-format   verify that the code formatter has been run"
	@echo "check-slow     run tests and validations that take awhile"
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
	poetry run python -m black --check -t py38 -l 80 polyA/ tests/

.PHONY: check-slow
check-slow:
	poetry run python -m pytest -m slow tests/ polyA/

.PHONY: format
format:
	poetry run python -m black -t py38 -l 80 polyA/ tests/

.PHONY: setup
setup:
	poetry install --no-dev

.PHONY: setup-dev
setup-dev:
	poetry install
