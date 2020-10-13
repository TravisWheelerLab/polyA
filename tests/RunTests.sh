#!/usr/bin/env bash

set -e

for f in ../fixtures/*.align.format ;
do echo "$f";
  python ../test_inputs/AdjudicateRegions.py --seqpos --lambda .1227 "$f" ../test_inputs/25p41g_edited.matrix;
  printf -- "-------------------------------------------------------------\n";
done
