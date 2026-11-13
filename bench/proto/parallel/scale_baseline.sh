#!/bin/bash
# baseline scaling: interleaved core counts per rep, so drift in background load hits all core counts
cd /Users/rvialle/suez
for rep in 1 2 3; do
  for ds in S M; do for c in 1 2 4 8; do bench/bench.sh baseline $ds $c p$rep; done; done
  for c in 1 4 8; do bench/bench.sh baseline M_dense $c p$rep; done
done
