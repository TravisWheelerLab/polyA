.PHONY: check-fast
check-fast:
	python -m pytest

.PHONY: check-format
check-format:
	python -m black --check .

.PHONY: format
format:
	python -m black .
