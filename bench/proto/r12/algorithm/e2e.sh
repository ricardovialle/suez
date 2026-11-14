#!/bin/bash
# interleaved end-to-end base12 vs algorithm_zz, 1 core
cd /Users/rvialle/suez
export RUNS=bench/proto/r12/algorithm/runs.tsv
for rep in 1 2 3; do for ds in $@; do
  bench/bench.sh algorithm_base $ds 1 $rep > /dev/null; bench/bench.sh algorithm_zz $ds 1 $rep > /dev/null
done; done
