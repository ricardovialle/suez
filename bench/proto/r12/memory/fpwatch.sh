#!/bin/bash
# Usage: fpwatch.sh <command...>: like bench/proto/review/memwatch.sh, but also samples (every ~0.3 s) the combined
# phys_footprint of the process tree with footprint(1) (private dirty + compressed; shared COW pages counted once).
"$@" & pid=$!
tree() { echo $1; for c in $(pgrep -P $1); do tree $c; done; }
peak=0; fpk=0; np=0
while kill -0 $pid 2>/dev/null; do
  pids=$(tree $pid)
  s=$(ps -o rss= -p $(echo $pids | tr ' ' ,) 2>/dev/null | awk '{s += $1} END {print s + 0}')
  [ "$s" -gt "$peak" ] && peak=$s
  args=""; n=0; for p in $pids; do args="$args -p $p"; n=$((n+1)); done
  f=$(footprint $args --noCategories -f bytes 2>/dev/null | awk '/Footprint:/ {for (i=1;i<=NF;i++) if ($i=="Footprint:") v=$(i+1)} END {print v+0}')
  # with several processes footprint ends with the combined total line; keep the last value printed
  [ "${f%.*}" -gt "$fpk" ] && { fpk=${f%.*}; np=$n; }
  sleep 0.2
done
wait $pid; rc=$?
echo "FPWATCH peak_tree_rss_mb=$((peak / 1024)) peak_tree_footprint_mb=$((fpk / 1048576)) procs_at_fp_peak=$np" >&2
exit $rc
