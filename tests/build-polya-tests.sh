#!/bin/sh

set -e

echo "tests:"

for f in fixtures/*.align.format; do
  echo "    $f:"
  echo "        command: python -m polyA --seqpos --lambda 0.1227 '$f' fixtures/25p41g_edited.matrix"
  echo "        stdout:"
  echo "            file: $f.output"
done
