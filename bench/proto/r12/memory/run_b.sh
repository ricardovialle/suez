#!/bin/bash
# interleaved: per-gene R overhead (P_many, S at 1 core) and multi-core RAM (M 2/4/8, XL_grm 4/8)
cd /Users/rvialle/suez; export RUNS=bench/proto/r12/memory/runs_b.tsv OUTDIR=bench/proto/r12/memory/out
B=bench/proto/r12/memory/bench_mem.sh
for rep in 1 2 3; do
  for ds in P_many S; do for impl in base12 memory_rside memory_all; do $B $impl $ds 1 $rep > /dev/null || echo FAIL $impl $ds 1 $rep; done; done
  for c in 2 4 8; do for impl in base12 memory_all; do $B $impl M $c $rep > /dev/null || echo FAIL $impl M $c $rep; done; done
  for c in 4 8; do for impl in base12 memory_all; do $B $impl XL_grm $c $rep > /dev/null || echo FAIL $impl XL_grm $c $rep; done; done
done
