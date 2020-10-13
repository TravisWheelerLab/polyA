#!/usr/bin/env bash

set -e

for f in ../fixtures/*.align.format ;
do echo "$f";
  python -m polyA --seqpos --lambda .1227 "$f" ../fixtures/25p41g_edited.matrix;
  printf -- "-------------------------------------------------------------\n";
done

#tests esl_scorematrix
for f in ../fixtures/*.align.format ;
do echo "$f";
  python -m polyA --seqpos "$f" ../fixtures/25p41g_edited.matrix;
  printf -- "-------------------------------------------------------------\n";
done

#tests soda output
for f in ../fixtures/*.align.format ;
do echo "$f";
  python -m polyA --seqpos --viz outfile.viz --heatmap outfile.heatmap "$f" ../fixtures/25p41g_edited.matrix;
  rm outfile.viz;
  rm outfile.viz.json;
  rm outfile.heatmap;
  printf -- "-------------------------------------------------------------\n";
done

# TODO: this is broken
##tests prior counts
#for f in ../fixtures/*.align.format ;
#do echo "$f";
#  python -m polyA --seqpos --lambda .1227 --priorCounts ../fixtures/SubfamPriorCounts.txt "$f" ../fixtures/25p41g_edited.matrix;
#  printf -- "-------------------------------------------------------------\n";
#done
