#!/bin/bash
# validate r9_rot against the r8_fix reference (1 core, rep 0)
cd /Users/rvialle/suez
for d in "$@"; do
  RUNS=bench/proto/algorithm/r9/runs_validate.tsv bench/bench.sh ${IMPL:-r9_rot} $d 1 0 > /dev/null || { echo "$d RUN FAILED"; continue; }
  echo -n "$d: "; Rscript bench/proto/algorithm/r9/cmp.R bench/results/out/r8_fix_${d}_step2_c1_r0.rds bench/results/out/${IMPL:-r9_rot}_${d}_step2_c1_r0.rds
done
