#!/usr/bin/env bash

set -e

for f in ../fixtures/ultra_test_files/*.fa ;
do echo "$f";
  ultra_output="${f}.ultra";
  align_file="${f}.cm.sto";
  m="${f}.cm.matrix";

  set -x
  python -m polyA --sequence-position --ultra-data "$ultra_output" "$align_file" "$m";
  set +x

  printf -- "-------------------------------------------------------------\n";
done


for f in ../fixtures/ultra_test_files/*.fa ;
do echo "$f";
  ultra_output="${f}.ultra";
  align_file="${f}.cm.sto";
  m="${f}.cm.matrix";

  set -x
  python -m polyA --sequence-position --soda --ultra-data "$ultra_output" "$align_file" "$m";
  set +x

  ls output.*
  rm output.*.viz;
  rm output.*.viz.json;
  printf -- "-------------------------------------------------------------\n";
done


for f in ../fixtures/ultra_test_files/*.fa ;
do echo "$f";
  ultra_output="${f}.ultra";
  align_file="${f}.cm.sto";
  m="${f}.cm.matrix";

  set -x
  python -m polyA --shard-gap 50 --sequence-position --ultra-data "$ultra_output" "$align_file" "$m";
  set +x

  printf -- "-------------------------------------------------------------\n";
done


for f in ../fixtures/ultra_test_files/*.fa ;
do echo "$f";
  ultra_output="${f}.ultra";
  align_file="${f}.cm.sto";
  m="${f}.cm.matrix";

  set -x
  python -m polyA --sequence-position --complexity-adjustment --ultra-data "$ultra_output" "$align_file" "$m";
  set +x

  printf -- "-------------------------------------------------------------\n";
done