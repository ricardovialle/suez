#!/bin/bash
# interleaved comparison of baseline and the prototypes, 3 reps
cd /Users/rvialle/suez; R=bench/proto/parallel/run_variant.sh
for rep in 1 2 3; do
  for ds in S M M_dense; do
    for c in 4 8; do for impl in baseline par_static par_dyn par_chunk4 par_chunk8; do $R $impl $ds $c v$rep; done; done
    for t in 4 8; do $R par_thr$t $ds 1 v$rep; done
    $R par_thr2 $ds 4 v$rep     # hybrid: 4 forked workers x 2 threads
  done
  for ds in S M; do $R baseline $ds 1 v$rep; $R par_thr1 $ds 1 v$rep; $R par_thr2 $ds 1 v$rep; done
done
