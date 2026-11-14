#!/bin/bash
# interleaved base12 vs memory_eig, 1 core
cd /Users/rvialle/suez; export RUNS=bench/proto/r12/memory/runs_eig.tsv
for rep in r1 r2 r3; do for ds in M XL XL_grm XXL; do for impl in base12 memory_eig; do
  OUTDIR=bench/proto/r12/memory/out bench/proto/r12/memory/bench_mem.sh $impl $ds 1 $rep > /dev/null || echo FAIL $impl $ds $rep
done; done; done
