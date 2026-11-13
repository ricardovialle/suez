#!/bin/bash
# Usage: bench/proto/memory/treerss.sh <command...>
# Runs the command and samples (every ~50-100 ms) its whole process tree (doMC forks included, which
# /usr/bin/time -l does not count). Prints to stderr (TREERSS line):
#   sum_mb      peak of the sum of RSS over the tree (copy-on-write pages shared by forks count once per
#               process: an upper bound),
#   max_proc_mb peak RSS of the largest single process,
#   fp_mb       peak of the tree's footprint (`footprint` Summary: private dirty + compressed memory,
#               shared pages counted once: the physical memory the run really costs),
#   nproc       peak number of processes.
"$@" &
root=$!
peak_sum=0; peak_max=0; peak_n=0; peak_fp=0
while kill -0 $root 2>/dev/null; do
  out=$(ps -A -o pid=,ppid=,rss= | awk -v r=$root '
    { pp[$1]=$2; rss[$1]=$3 }
    END { in_[r]=1; changed=1
          while (changed) { changed=0; for (p in pp) if (!(p in in_) && (pp[p] in in_)) { in_[p]=1; changed=1 } }
          s=0; m=0; n=0; l=""; for (p in in_) if (p in rss) { s+=rss[p]; n++; if (rss[p]>m) m=rss[p]; l=l" -p "p }
          print s, m, n, l }')
  read s m n pids <<< "$out"
  [ "${s:-0}" -gt $peak_sum ] && peak_sum=$s
  [ "${m:-0}" -gt $peak_max ] && peak_max=$m
  [ "${n:-0}" -gt $peak_n ] && peak_n=$n
  f=$(footprint $pids 2>/dev/null | awk '
    /Summary Footprint:/ { v=$3; u=$4 } /^[^ ].*Footprint:/ && !/Summary/ && v=="" { for (i=1;i<=NF;i++) if ($i=="Footprint:") { v1=$(i+1); u1=$(i+2) } }
    END { if (v=="") { v=v1; u=u1 } if (u=="MB") v*=1024; if (u=="GB") v*=1048576; printf "%d", v }')
  [ "${f:-0}" -gt $peak_fp ] && peak_fp=$f
  sleep 0.03
done
wait $root; rc=$?
awk -v s=$peak_sum -v m=$peak_max -v n=$peak_n -v f=$peak_fp 'BEGIN { printf "TREERSS\tsum_mb=%.1f\tmax_proc_mb=%.1f\tfp_mb=%.1f\tnproc=%d\n", s/1024, m/1024, f/1024, n }' >&2
exit $rc
