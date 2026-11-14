#!/bin/bash
# Round 12+ benchmark: implementations interleaved within each repetition, serial.
#   bench/round12_bench.sh <reps> <impl> [impl ...]     (impls = builds in bench/lib/<impl>)
# Env: RUNS (output tsv, default bench/results/runs_r12.tsv), DS1 (1-core datasets), DSC (multi-core datasets),
#      CORES (multi-core counts), STEP1 (datasets for learn_covariance).
reps=$1; shift; impls_fwd="$*"
export RUNS=${RUNS:-bench/results/runs_r12.tsv}
DS1=${DS1:-"S M L XL XXL M_grm L_grm XL_grm M_learn M_lc M_num M_fac3 M_eqtl M_cov M_cov30 M_miss M_dense M_unbal S_rare HS LA P_many P_many10 R_real XL_grm_many"}
DSC=${DSC:-"M M_dense XL_grm R_real"}
CORES=${CORES:-"2 4 8"}
STEP1=${STEP1:-"S M"}
for r in $(seq 1 $reps); do
  # alternate the implementation order between repetitions (position bias, reviewer A5)
  if [ $((r % 2)) -eq 0 ]; then impls=$(echo $impls_fwd | tr ' ' '\n' | tail -r | tr '\n' ' '); else impls=$impls_fwd; fi
  for d in $DS1; do for i in $impls; do bench/bench.sh $i $d 1 r$r > /dev/null || echo "FAILED $i $d 1 r$r"; done; done
  for d in $DSC; do for c in $CORES; do for i in $impls; do bench/bench.sh $i $d $c r$r > /dev/null || echo "FAILED $i $d $c r$r"; done; done; done
  for d in $STEP1; do for i in $impls; do bench/bench.sh $i $d 1 r$r step1 > /dev/null || echo "FAILED $i $d step1 r$r"; done; done
done
