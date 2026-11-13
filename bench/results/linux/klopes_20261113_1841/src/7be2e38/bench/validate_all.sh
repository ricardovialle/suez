#!/bin/bash
# End-to-end statistical validation of a build against the baseline on every dataset (1 core, rep 0).
# Usage: bench/validate_all.sh <impl>  -> bench/results/validate_<impl>.txt
impl=$1; out=bench/results/validate_$impl.txt; : > $out
export RUNS=bench/results/runs_validate.tsv
for d in S M L XL M_grm M_learn M_num M_fac3 M_eqtl M_cov M_miss M_dense S_rare LA L_grm M_unbal M_lc HS; do
  [ -f bench/results/out/baseline_${d}_step2_c1_r0.rds ] || bench/bench.sh baseline $d 1 0 > /dev/null
  bench/bench.sh $impl $d 1 0 > /dev/null
  echo "=== $d" >> $out
  Rscript bench/proto/review/validate.R compare bench/results/out/baseline_${d}_step2_c1_r0.rds bench/results/out/${impl}_${d}_step2_c1_r0.rds $d >> $out 2>&1
done
grep -E '^===|VERDICT' $out
