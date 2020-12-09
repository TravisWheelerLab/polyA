#!/usr/bin/env bash

set -e

# test with ultra output
for f in ../fixtures/ultra_test_files/*.fa ;
do echo "$f";
  ultra_output="${f}.ultra";
  align_file="${f}.cm.sto";
  m="${f}.cm.matrix";
  python -m polyA --seq-pos --ultra-output "$ultra_output" "$align_file" "$m";
  printf -- "-------------------------------------------------------------\n";
done


#tests soda output
for f in ../fixtures/ultra_test_files/*.fa ;
do echo "$f";
  ultra_output="${f}.ultra";
  align_file="${f}.cm.sto";
   m="${f}.cm.matrix";
  python -m polyA --seq-pos --soda --heatmap --ultra-output "$ultra_output" "$align_file" "$m"";
  ls polya-output.*
  rm polya-output.*.viz;
  rm polya-output.*.viz.json;
  rm polya-output.*.heatmap;
  printf -- "-------------------------------------------------------------\n";
done