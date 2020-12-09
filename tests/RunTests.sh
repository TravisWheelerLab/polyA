#!/usr/bin/env bash

set -e

for f in ../fixtures/ex*.sto;
do echo "$f";
  m=${f%sto}matrix
  python -m polyA --sequence-position "$f" "$m";
  printf -- "-------------------------------------------------------------\n";
done

## tests confidence only option - use --lambda just so doesn't run esl_scorematrix
#for f in ../fixtures/ex*.cm.sto;
#do echo "$f";
#  m=${f%sto}matrix
#  python -m polyA --confidence "$f" "$m";
#  printf -- "-------------------------------------------------------------\n";
#done
#
##tests soda output
#for f in ../fixtures/ex*.cm.sto;
#do echo "$f";
#  m=${f%sto}matrix
#  python -m polyA --sequence-position --soda --heatmap "$f" "$m";
#  ls output.*
#  rm output.*.viz;
#  rm output.*.viz.json;
#  rm output.*.heatmap;
#  printf -- "-------------------------------------------------------------\n";
#done
#
###tests prior counts
#for f in ../fixtures/ex*.cm.sto;
#do echo "$f";
#  m=${f%sto}matrix
#  python -m polyA --sequence-position --prior-counts ../fixtures/SubfamPriorCounts.txt "$f" "$m";
#  printf -- "-------------------------------------------------------------\n";
#done
