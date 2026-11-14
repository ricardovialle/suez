#!/bin/bash
# interleaved large-N 1 core: base12 vs memory_all2 (lean set-up + R side) vs memory_exact (BLAS-exact set-up + R side)
cd /Users/rvialle/suez; export RUNS=bench/proto/r12/memory/runs_c.tsv OUTDIR=bench/proto/r12/memory/out
for rep in 1 2 3; do for ds in XL XL_grm XXL; do for impl in base12 memory_all2 memory_exact; do
  bench/proto/r12/memory/bench_mem.sh $impl $ds 1 c$rep > /dev/null || echo FAIL $impl $ds $rep
done; done; done
