#!/bin/bash
# Kernel-level timing, process-interleaved: for each rep, each dataset, base12 / integ_c1 / integ_c2 in turn.
#   bench/proto/r12/numerical/integ_timing.sh <reps> <out.tsv> ds:genes ...
reps=$1; out=$2; shift 2; SP=$(mktemp -d)
echo -e "rep\tdataset\tlib\ttime" > $out
for r in $(seq $reps); do for spec in "$@"; do ds=${spec%%:*}; ng=${spec##*:}
  for L in ${LIBS:-base12 integ_c1 integ_c2}; do
    t=$(Rscript bench/proto/r12/numerical/integ_run.R bench/lib/$L $ds $ng $SP/${L}_${ds}.rds 2>/dev/null | grep ' time ' | awk '{print $NF}')
    echo -e "$r\t$ds\t$L\t$t" >> $out
  done
  if false; then :
    SUEZ_HARNESS_LIB=bench/lib/base12 Rscript bench/proto/r12/numerical/cmp_rds.R $SP/integ_c1_${ds}.rds $SP/integ_c2_${ds}.rds 2>/dev/null | grep identical | sed "s/^/$ds c2 vs c1: /"; fi
done; done
