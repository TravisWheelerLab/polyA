#!/usr/bin/env bash

set -e

for f in ../fixtures/*.align.format ;
do echo "$f";
  python -m polyA --seqpos --lambda .1227 "$f" ../fixtures/25p41g_edited.matrix;
  printf -- "-------------------------------------------------------------\n";
done
