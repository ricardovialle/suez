#!/bin/bash
# interleaved base12 / integ_c2 / integ_c3, 1 core, 3 reps; then the reviewer's mem_gene.R (N = 2000, S = 0 / 4000)
cd /Users/rvialle/suez; export RUNS=bench/proto/r12/memory/runs_c3.tsv OUTDIR=bench/proto/r12/memory/out
for rep in 1 2 3; do for ds in M XL_grm XL_grm_many XXL R_real; do for impl in base12 integ_c2 integ_c3; do
  bench/proto/r12/memory/bench_mem.sh $impl $ds 1 g$rep > /dev/null || echo FAIL $impl $ds $rep
done; done; done
for rep in 1 2 3; do for S in 0 4000; do for l in base12 integ_c2 integ_c3; do
  /usr/bin/time -l Rscript bench/proto/r12/review/mem_gene.R bench/lib/$l 1000 $S 2>&1 | awk -v l=$l -v S=$S -v r=$rep '/fit/ {f=$5} /maximum resident/ {m=$1/1048576} /peak memory footprint/ {p=$1/1048576} END {printf "%s\t%s\t%s\t%.0f\t%.0f\t%s\n", l, S, r, m, p, f}' >> bench/proto/r12/memory/mem_gene_c3.tsv
done; done; done
