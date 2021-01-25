#!/usr/bin/env bash

set -e

for f in ../fixtures/ex*.sto;
do
  m=${f%sto}matrix

  set -x
  python -m polyA "$f" "$m";
  set +x

  printf -- "-------------------------------------------------------------\n";
done

for f in ../fixtures/ex*.sto;
do
  m=${f%sto}matrix

  set -x
  python -m polyA --confidence "$f" "$m";
  set +x

  printf -- "-------------------------------------------------------------\n";
done

for f in ../fixtures/ex*.sto;
do
  m=${f%sto}matrix

  set -x
  python -m polyA --sequence-position --soda --heatmap "$f" "$m";
  set +x

  ls output.*
  rm output.*.viz;
  rm output.*.viz.json;
  rm output.*.heatmap;
  printf -- "-------------------------------------------------------------\n";
done

for f in ../fixtures/ex*cm.sto;
do
  m=${f%sto}matrix

  set -x
  python -m polyA --sequence-position --prior-counts ../fixtures/SubfamPriorCounts.txt "$f" "$m";
  set +x

  printf -- "-------------------------------------------------------------\n";
done

for f in ../fixtures/ex*.sto;
do
  m=${f%sto}matrix

  set -x
  python -m polyA --matrix-position "$f" "$m";
  set +x

  printf -- "-------------------------------------------------------------\n";
done
