#!/usr/bin/env bash

set -e

for f in ../fixtures/ex*.cm.sto;
do echo "$f";
  m=${f%sto}matrix
  python -m polyA --seq-pos "$f" "$m";
  printf -- "-------------------------------------------------------------\n";
done

# tests confidence only option - use --lambda just so doesn't run esl_scorematrix
for f in ../fixtures/ex*.cm.sto;
do echo "$f";
  m=${f%sto}matrix
  python -m polyA --confidence "$f" "$m";
  printf -- "-------------------------------------------------------------\n";
done

#TODO: these will fail on git because no esl_scorematrix
##tests esl_scorematrix
#for f in ../fixtures/*.align.format ;
#do echo "$f";
#  python -m polyA --seqpos "$f" ../fixtures/25p41g_edited.matrix;
#  printf -- "-------------------------------------------------------------\n";
#done

#tests soda output
for f in ../fixtures/ex*.cm.sto;
do echo "$f";
  m=${f%sto}matrix
  python -m polyA --seq-pos --soda --heatmap "$f" "$m";
  ls polya-output.*
  rm polya-output.*.viz;
  rm polya-output.*.viz.json;
  rm polya-output.*.heatmap;
  printf -- "-------------------------------------------------------------\n";
done

##tests prior counts
for f in ../fixtures/ex*.cm.sto;
do echo "$f";
  m=${f%sto}matrix
  python -m polyA --seq-pos --prior-counts ../fixtures/SubfamPriorCounts.txt "$f" "$m";
  printf -- "-------------------------------------------------------------\n";
done
