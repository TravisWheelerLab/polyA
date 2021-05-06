#!/usr/bin/env bash
# given an alignment file and alignment tool, call the correct stockholm parser

if [ $# -lt 2 ]
then
  echo "$0: takes 2 arguments: alignment file, alignment tool" 1>&2
  exit 1
fi

alignment_file="$1"
alignment_tool="$2"

if [ "$alignment_tool" = "cross_match" ]; then
  python cm_to_stockholm.py "$alignment_file"
fi