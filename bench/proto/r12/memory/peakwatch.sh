#!/bin/bash
# peakwatch.sh <command...>: peak RSS of the command (sampled every 20 ms), printed as PEAK <MB>
"$@" & pid=$!; pk=0
while kill -0 $pid 2>/dev/null; do s=$(ps -o rss= -p $pid 2>/dev/null | tr -d ' '); [ -n "$s" ] && [ "$s" -gt "$pk" ] && pk=$s; sleep 0.02; done
wait $pid; echo "PEAK $((pk / 1024))"
