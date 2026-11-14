#!/bin/bash
# Usage: bench/bench.sh <impl> <dataset> <cores> <rep> [step]
# Runs bench/run.R with the suez build in bench/lib/<impl> under /usr/bin/time -l, appends one row to
# bench/results/runs.tsv: impl dataset step cores rep wall_s user_s sys_s cpu_util n_tests maxrss_mb
# peak_footprint_mb blocks_in blocks_out tree_rss_mb. user/sys are from /usr/bin/time (forked workers included,
# R start-up included); tree_rss_mb = peak RSS summed over the process tree (bench/proto/review/memwatch.sh,
# an upper bound: shared pages count once per fork); and keeps the result table in bench/results/out/.
set -e
impl=$1; ds=$2; cores=$3; rep=$4; step=${5:-step2}
mkdir -p ${OUTDIR:-bench/results/out}
out=${OUTDIR:-bench/results/out}/${impl}_${ds}_${step}_c${cores}_r${rep}.rds
log=$(mktemp)
OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=$cores \
  bench/proto/review/memwatch.sh /usr/bin/time -l Rscript bench/run.R bench/lib/$impl $ds $cores $out $step > $log.out 2> $log
r=$(grep '^RESULT' $log.out | cut -f2-)
rss=$(awk '/maximum resident set size/ {print $1/1048576}' $log)
pf=$(awk '/peak memory footprint/ {print $1/1048576}' $log)
bi=$(awk '/block input operations/ {print $1}' $log)
bo=$(awk '/block output operations/ {print $1}' $log)
tu=$(awk '/ real .* user .* sys/ {print $3}' $log); ts=$(awk '/ real .* user .* sys/ {print $5}' $log)
tree=$(sed -n 's/.*peak_tree_rss_mb=\([0-9]*\).*/\1/p' $log)
[ -z "$r" ] && { cat $log $log.out | tail -20; exit 1; }
f=${RUNS:-bench/results/runs.tsv}
[ -f $f ] || printf 'impl\tdataset\tstep\tcores\trep\twall_s\tuser_s\tsys_s\tcpu_util\tn_tests\tmaxrss_mb\tpeak_footprint_mb\tblocks_in\tblocks_out\ttree_rss_mb\n' > $f
set -- $r
util=$(awk -v w=$1 -v u=$tu -v s=$ts 'BEGIN{printf "%.2f", (u+s)/w}')
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.1f\t%.1f\t%s\t%s\t%s\n' $impl $ds $step $cores $rep $1 $tu $ts $util $4 $rss $pf $bi $bo $tree | tee -a $f
rm -f $log $log.out
