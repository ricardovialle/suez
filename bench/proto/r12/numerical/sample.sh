#!/bin/bash
# sample.sh <proto.cpp> <dataset> <genes> <out.txt>: macOS `sample` of the prototype's fitting loop (after set-up)
SP=$(dirname $4); Rscript bench/proto/r12/numerical/loop.R $1 $2 $3 200 > $4.log 2>&1 &
until grep -q READY $4.log; do sleep 1; done; sleep 2
PID=$(pgrep -f "loop.R --args $1 $2" | head -1); sample $PID 8 -file $4 >/dev/null 2>&1; kill $PID
grep -A30 'Sort by top of stack' $4 | tail -29 | grep -v -E '__workq|__semwait|^ *$' | sed -E 's/\(in sourceCpp[^)]*\)//; s/Eigen::internal:://g; s/Eigen::Matrix<double, -1, -1, 0, -1, -1>/MXd/g; s/  +/ /g' | awk '{n=$NF; $NF=""; print n, substr($0,1,100)}' | head -22
