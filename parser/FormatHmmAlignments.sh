#!/usr/bin/env bash

for f in ../fixtures/hmm_test_files/*.fa.hmm_align ;
do echo "$f";
  python hmm_align_to_table.py $f;
  table_output="${f}.tbl.txt";
  sort -k6 -n $table_output > $f;
  rm $table_output
  python hmm_table_to_stockholm.py $f;
done