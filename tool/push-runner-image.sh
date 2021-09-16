#!/usr/bin/env sh

set -e

POLYA_VERSION=$(python polyA/_version.py)

docker push traviswheelerlab/polya:${POLYA_VERSION}
docker push traviswheelerlab/polya:latest
