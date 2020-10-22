#!/usr/bin/env bash

set -e

# test with ultra output
for f in ../fixtures/ultra_test_files/*.fa ;
do echo "$f";
  ultra_output="${f}.ultra";
  align_file="${f}.cm.align";
  python -m polyA --seqpos --lambda .1227 --ultraOutput "$ultra_output" "$align_file" ../fixtures/25p41g_edited.matrix;
  printf -- "-------------------------------------------------------------\n";
done

#tests soda output
for f in ../fixtures/ultra_test_files/*.fa ;
do echo "$f";
  ultra_output="${f}.ultra";
  align_file="${f}.cm.align";
  python -m polyA --seqpos --lambda .1227 --viz outfile.viz --heatmap outfile.heatmap --ultraOutput "$ultra_output" "$align_file" ../fixtures/25p41g_edited.matrix;
  rm outfile.viz;
  rm outfile.viz.json;
  rm outfile.heatmap;
  printf -- "-------------------------------------------------------------\n";
done
