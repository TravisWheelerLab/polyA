#!/usr/bin/env bash

set -e

for f in ../fixtures/ex*.sto;
do
  m=${f%sto}matrix

  set -x
  python -m polyA --benchmark "$f.report" "$f" "$m"
  set +x

  printf -- "-------------------------------------------------------------\n";
done
