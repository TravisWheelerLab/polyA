#!/usr/bin/env sh

set -e

BENCHMARK_BASE=benchmarks
BENCHMARK_FILE=chr8_debug_large

BENCHMARK_INPUT=$BENCHMARK_BASE/$BENCHMARK_FILE

pipenv run python -u -m scalene \
  --outfile $BENCHMARK_INPUT.profile.html \
  --html \
  --program-path . \
  --cpu-only \
  --reduced-profile \
  $BENCHMARK_BASE/entrypoint.py \
    --shard-gap 0 \
    --log-level verbose \
    --sequences $BENCHMARK_INPUT.fa \
    $BENCHMARK_INPUT.fa.cm.sto \
    $BENCHMARK_INPUT.fa.cm.matrix \
| tee $BENCHMARK_INPUT.output

# Determine whether the output has changed from our gold file.
diff -q \
  $BENCHMARK_INPUT.gold \
  $BENCHMARK_INPUT.output
