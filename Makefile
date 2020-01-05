.PHONY: check
check: check-fast check-slow check-format

.PHONY: check-fast
check-fast:
	pipenv run python -m pytest -m 'not slow'

.PHONY: check-format
check-format:
	pipenv run python -m black --check polyA/ tests/

.PHONY: check-slow
check-slow:
	pipenv run python -m pytest -m slow

.PHONY: format
format:
	pipenv run python -m black polyA/ tests/

.PHONY: setup
setup:
	pipenv install

.PHONY: setup-dev
setup-dev:
	pipenv install --dev
