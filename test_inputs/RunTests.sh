#!/usr/bin/env bash

## ---------
## UTILITIES
## ---------

RED='\033[0;31m'
GREEN='\033[0;32m'
NOCOLOR='\033[00m'

BALLOTX='✘'
CHECKMARK='✔'

function assert-same() {
  local base_file="$1"
  local pl_file="$base_file.pl.output"
  local py_file="$base_file.py.output"
  local diff_file="$base_file.output.diff"
  local output_diff
  output_diff=$(diff "$pl_file" "$py_file")
  if [ "$output_diff" == "" ]; then
    echo -e "${GREEN}${CHECKMARK} $base_file ... PASSED${NOCOLOR}"
  else
    echo -e "${RED}${BALLOTX} $base_file ... FAILED${NOCOLOR} (see $diff_file)"
    echo "$output_diff" > "$diff_file"
  fi
  rm "$pl_file"
  rm "$py_file"
}

function run-script() {
  local alignment="$1"
  local interpreter="$2"
  local extension="$3"
  "$interpreter" "AdjudicateRegions.$extension" \
      --lambda 0.1227 \
      "$alignment" \
      25p41g_edited.matrix \
      > "$alignment.$extension.output" \
      2>&1
}

function run-perl() {
  run-script "$1" "perl" "pl"
}

function run-python() {
  run-script "$1" "python" "py"
}

## -----
## TESTS
## -----

# artificial test files most often used
for alignment in *.align; do
  run-perl "$alignment"
  run-python "$alignment"
  assert-same "$alignment"
done

# real test files from ISB, not as frequently used but need to use later 
for alignment in *.align.format; do
  run-perl "$alignment"
  run-python "$alignment"
  assert-same "$alignment"
done
>>>>>>> 3bd684e84c74e8f915cf824ad61cc6873b327869
