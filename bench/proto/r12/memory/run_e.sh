#!/bin/bash
# interleaved: memory_gene (all memory changes) vs base12; explicit-kernel large N, P_many10, M; XL_grm 8 cores
cd /Users/rvialle/suez; export RUNS=bench/proto/r12/memory/runs_e.tsv OUTDIR=bench/proto/r12/memory/out
B=bench/proto/r12/memory/bench_mem.sh
for rep in 1 2 3; do
  for ds in M P_many10 XL_grm XL_grm_many; do for impl in base12 memory_gene; do $B $impl $ds 1 e$rep > /dev/null || echo FAIL $impl $ds; done; done
  for impl in base12 memory_gene; do $B $impl XL_grm 8 e$rep > /dev/null || echo FAIL $impl XL_grm 8; done
done
