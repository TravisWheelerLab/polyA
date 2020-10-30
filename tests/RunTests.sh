#!/usr/bin/env bash

set -e

for f in ../fixtures/ex*.sto;
do echo "$f";
  python -m polyA --seqpos --lambda .1227 "$f" ../fixtures/25p41g_edited.matrix;
  printf -- "-------------------------------------------------------------\n";
done

# tests confidence only option
for f in ../fixtures/ex*.sto;
do echo "$f";
  python -m polyA --confidence "$f" ../fixtures/25p41g_edited.matrix;
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
for f in ../fixtures/ex*.sto;
do echo "$f";
  python -m polyA --seqpos --lambda .1227 --viz outfile.viz --heatmap outfile.heatmap "$f" ../fixtures/25p41g_edited.matrix;
  rm outfile.viz;
  rm outfile.viz.json;
  rm outfile.heatmap;
  printf -- "-------------------------------------------------------------\n";
done

##tests prior counts
for f in ../fixtures/ex*.sto;
do echo "$f";
  python -m polyA --seqpos --lambda .1227 --priorCounts ../fixtures/SubfamPriorCounts.txt "$f" ../fixtures/25p41g_edited.matrix;
  printf -- "-------------------------------------------------------------\n";
done
