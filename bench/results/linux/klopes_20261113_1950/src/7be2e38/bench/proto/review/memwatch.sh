#!/bin/bash
# Usage: bench/proto/review/memwatch.sh <command...>
# Runs the command and samples (every 0.1 s) the summed RSS of its whole process tree (doMC/mclapply forks
# included), which /usr/bin/time -l does not report: its maxrss is the max of single processes, its peak
# footprint is the direct child only. Summed RSS double-counts copy-on-write pages shared by the forks
# (~220 MB of loaded libraries each), so it is an upper bound; report it next to /usr/bin/time's maxrss.
# shortcut: polling misses peaks shorter than 0.1 s; fine for the multi-second benchmark runs.
"$@" & pid=$!
tree() { echo $1; for c in $(pgrep -P $1); do tree $c; done; }
peak=0; np=0
while kill -0 $pid 2>/dev/null; do
  pids=$(tree $pid | paste -sd, -)
  s=$(ps -o rss= -p $pids 2>/dev/null | awk '{s += $1} END {print s + 0}')
  if [ "$s" -gt "$peak" ]; then peak=$s; np=$(echo $pids | tr , '\n' | wc -l | tr -d ' '); fi
  sleep 0.1
done
wait $pid; rc=$?
echo "MEMWATCH peak_tree_rss_mb=$((peak / 1024)) procs_at_peak=$np" >&2
exit $rc
