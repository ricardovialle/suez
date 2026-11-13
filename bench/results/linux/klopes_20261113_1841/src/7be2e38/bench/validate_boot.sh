#!/bin/bash
# validate_all.sh plus the bootstrap distribution check: bench/validate_boot.sh <impl>
impl=$1; bench/validate_all.sh $impl > /dev/null
out=bench/results/validate_boot_$impl.txt; : > $out
for d in S M L XL M_grm M_learn M_num M_fac3 M_eqtl M_cov M_miss M_dense S_rare LA L_grm M_unbal M_lc HS; do
  printf '%s\t' $d >> $out
  Rscript bench/boot_distribution.R bench/results/out/baseline_${d}_step2_c1_r0.rds bench/results/out/${impl}_${d}_step2_c1_r0.rds >> $out 2>&1
done
