#!/bin/bash
# interleaved end-to-end timing r8_fix vs r9_rot
cd /Users/rvialle/suez; export RUNS=bench/proto/algorithm/r9/runs_timing${SUF}.tsv
for rep in ${REPS:-11 12 13}; do
  for d in M L XL M_dense; do for impl in ${BASE:-r8_fix} ${IMPL:-r9_rot}; do bench/bench.sh $impl $d 1 $rep > /dev/null; done; done
  for d in M M_dense; do for impl in ${BASE:-r8_fix} ${IMPL:-r9_rot}; do bench/bench.sh $impl $d 4 $rep > /dev/null; done; done
done
