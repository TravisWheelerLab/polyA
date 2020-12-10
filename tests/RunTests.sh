#!/usr/bin/env bash

set -e

python -m polyA ../fixtures/ex1.fa.align.sto ../fixtures/ex1.fa.align.matrix
python -m polyA ../fixtures/ex4.fa.align.sto ../fixtures/ex4.fa.align.matrix
