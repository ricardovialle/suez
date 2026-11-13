#!/bin/bash
# memsample.sh <impl> <ds> <cores>: runs bench/run.R and samples the process tree every 0.25 s; prints the peak
# number of R processes, peak summed RSS and peak summed phys_footprint (private dirty, via footprint(1)).
cd /Users/rvialle/suez; impl=$1
case $impl in par_thr*) export SUEZ_THREADS=${impl#par_thr};; par_chunk4) export SUEZ_SCHEDULE=chunk;; esac
OMP_NUM_THREADS=1 Rscript bench/run.R bench/lib/$impl $2 $3 /dev/null & root=$!
maxn=0; maxrss=0; maxfp=0
while kill -0 $root 2>/dev/null; do
  pids="$root $(pgrep -P $root | tr '\n' ' ')"
  rss=$(ps -o rss= -p ${pids// /,} 2>/dev/null | awk '{s+=$1} END {print s/1024}')
  fp=$(footprint $pids 2>/dev/null | awk '/phys_footprint:/ {v=$2; u=$3; if (u=="KB") v/=1024; if (u=="GB") v*=1024; s+=v} END {print s+0}')
  n=$(echo $pids | wc -w)
  maxn=$(( n > maxn ? n : maxn )); maxrss=$(awk -v a=$rss -v b=$maxrss 'BEGIN{print (a>b)?a:b}'); maxfp=$(awk -v a=$fp -v b=$maxfp 'BEGIN{print (a>b)?a:b}')
  sleep 0.25
done
echo -e "$impl\t$2\tcores=$3\tprocs=$maxn\tsum_rss_mb=$maxrss\tsum_footprint_mb=$maxfp"
