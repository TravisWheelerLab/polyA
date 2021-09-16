#!/usr/bin/env sh

set -e

POLYA_VERSION=$(python polyA/_version.py)

docker build -f Dockerfile_run \
    -t traviswheelerlab/polya:${POLYA_VERSION} \
    -t traviswheelerlab/polya:latest \
    $@ \
    .
