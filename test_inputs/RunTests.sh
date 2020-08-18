#!/usr/bin/env bash

set -e

for f in *.align.format ;
do echo "$f";
  poetry run python3 AdjudicateRegions.py --seqpos --lambda .1227 $f 25p41g_edited.matrix;
  printf -- "-------------------------------------------------------------\n";
done
