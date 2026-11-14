#!/bin/bash
# End-to-end runs (bench/bench.sh), interleaved: rep outer, configurations inner.
#   bench/proto/r12/parallel/e2e.sh <reps> <dataset> "<impl>:<cores>:<threads>:<sched>" ...
# impl = library under bench/lib (parallel_thr = threaded engine + scheduling; parallel_acc = the same + Accelerate
# eigen solver). Rows are labelled <impl>_t<threads>[_<sched>] through symlinks in bench/lib.
cd "$(dirname "$0")/../../../.." || exit 1
export RUNS=${RUNS:-bench/proto/r12/parallel/runs.tsv}
reps=$1; ds=$2; shift 2
for r in $(seq 1 $reps); do
  for c in "$@"; do
    IFS=: read impl cores thr sched <<< "$c"
    label=${impl}_t${thr}${sched:+_$sched}; [ "$impl" = base12 ] && label=parallel_base12   # own output names (bench/results/out)
    case $label in parallel_*) ;; *) label=parallel_$label ;; esac
    [ -e bench/lib/$label ] || ln -s $impl bench/lib/$label
    SUEZ_THREADS=$thr SUEZ_SCHED=${sched:-static} bench/bench.sh $label $ds $cores r$r > /dev/null || echo "FAILED $label $ds $cores"
  done
done
