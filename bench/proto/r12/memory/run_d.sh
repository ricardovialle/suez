#!/bin/bash
# interleaved: final variant memory_all3 (own bind, no dplyr) vs base12: P_many 1 core, M 2/4/8, XL_grm 8
cd /Users/rvialle/suez; export RUNS=bench/proto/r12/memory/runs_d.tsv OUTDIR=bench/proto/r12/memory/out
B=bench/proto/r12/memory/bench_mem.sh
for rep in 1 2 3; do
  for impl in base12 memory_all3; do $B $impl P_many 1 d$rep > /dev/null || echo FAIL; done
  for c in 2 4 8; do for impl in base12 memory_all3; do $B $impl M $c d$rep > /dev/null || echo FAIL; done; done
  for impl in base12 memory_all3; do $B $impl XL_grm 8 d$rep > /dev/null || echo FAIL; done
done
