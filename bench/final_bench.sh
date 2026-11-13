#!/bin/bash
# Final authoritative benchmark: baseline (tag perf-baseline) vs final (branch HEAD), serial on an idle machine,
# implementations interleaved within each repetition. Rows in bench/results/runs_final.tsv.
#   bench/final_bench.sh [reps]
reps=${1:-3}; export RUNS=bench/results/runs_final.tsv
DS="S M L XL M_grm M_learn M_num M_fac3 M_eqtl M_cov M_miss M_dense S_rare LA L_grm M_unbal M_lc HS"
for r in $(seq 1 $reps); do
  for d in $DS; do for i in baseline final; do bench/bench.sh $i $d 1 f$r > /dev/null; done; done
  for d in M M_dense; do for c in 2 4 8; do for i in baseline final; do bench/bench.sh $i $d $c f$r > /dev/null; done; done; done
  for d in S M; do for i in baseline final; do bench/bench.sh $i $d 1 f$r step1 > /dev/null; done; done
done
