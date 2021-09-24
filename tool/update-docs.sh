#!/usr/bin/env sh

set -e

sphinx-apidoc -f -o docs/ .
cd docs && make html
