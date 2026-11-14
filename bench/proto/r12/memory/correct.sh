#!/bin/bash
# outputs for correctness: Rscript bench/run.R per impl/dataset -> bench/proto/r12/memory/out/<impl>_<ds>_step2_c<cores>_ck.rds
cd /Users/rvialle/suez; cores=${CORES:-1}
for ds in $DSETS; do for impl in $IMPLS; do
  o=bench/proto/r12/memory/out/${impl}_${ds}_step2_c${cores}_ck.rds
  [ -f $o ] || OMP_NUM_THREADS=1 Rscript bench/run.R bench/lib/$impl $ds $cores $o > /dev/null 2>&1 || echo FAIL $impl $ds
done; done
