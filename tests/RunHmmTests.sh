#!/usr/bin/env bash

set -e

for f in ../fixtures/hmm_test_files/ex*.sto;
do echo "$f";
  set -x
  python -m polyA --sequence-position --hmm-path  ../fixtures/hmm_test_files/families.hmm.txt "$f" "";
  set +x

  printf -- "-------------------------------------------------------------\n";
done