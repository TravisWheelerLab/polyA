#!/usr/bin/env bash

set -e

for f in ../fixtures/*.align.format ;
do echo "$f";
  python AdjudicateRegions.py --seqpos --lambda .1227 "$f" 25p41g_edited.matrix;
  printf -- "-------------------------------------------------------------\n";
done
