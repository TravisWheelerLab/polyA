#!/usr/bin/env bash

set -e

# test with ultra output
for f in ../fixtures/ultra_test_files/*.fa ;
do echo "$f";
  ultra_output="${f}.ultra";
  align_file="${f}.cm.sto";
  python -m polyA --seq-pos --lambda .1227 --ultra-output "$ultra_output" "$align_file" ../fixtures/25p41g_edited.matrix;
  printf -- "-------------------------------------------------------------\n";
done

#tests soda output
for f in ../fixtures/ultra_test_files/*.fa ;
do echo "$f";
  ultra_output="${f}.ultra";
  align_file="${f}.cm.sto";
  python -m polyA --seq-pos --lambda .1227 --soda --heatmap --ultra-output "$ultra_output" "$align_file" ../fixtures/25p41g_edited.matrix;
  ls polya-output.*
  rm polya-output.*.viz;
  rm polya-output.*.viz.json;
  rm polya-output.*.heatmap;
  printf -- "-------------------------------------------------------------\n";
done