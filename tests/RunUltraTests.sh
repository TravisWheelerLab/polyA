#!/usr/bin/env bash

set -e

# test with ultra output
for f in ../fixtures/ultra_test_files/*.fa ;
do echo "$f";
  ultra_output="${f}.ultra";
  align_file="${f}.cm.sto";
  m="${f}.cm.matrix";
  python -m polyA --sequence-position --ultra-data "$ultra_output" "$align_file" "$m";
  printf -- "-------------------------------------------------------------\n";
done


#tests soda output
for f in ../fixtures/ultra_test_files/*.fa ;
do echo "$f";
  ultra_output="${f}.ultra";
  align_file="${f}.cm.sto";
   m="${f}.cm.matrix";
  python -m polyA --sequence-position --soda --heatmap --ultra-data "$ultra_output" "$align_file" "$m";
  ls output.*
  rm output.*.viz;
  rm output.*.viz.json;
  rm output.*.heatmap;
  printf -- "-------------------------------------------------------------\n";
done