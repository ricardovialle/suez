#!/bin/bash
# after validate12 finishes: tests/*.R against integ_c3, then the c3 benchmarks
cd /Users/rvialle/suez; S=/private/tmp/claude-501/-Users-rvialle-suez/73f927f3-8032-4547-afc6-66e11ee5023b/scratchpad
until grep -q "VERDICT" bench/proto/r12/memory/validate12_c3_vs_c2.log 2>/dev/null; do sleep 20; done
for d in tests/*.R; do Rscript -e ".libPaths(c('bench/lib/integ_c3', .libPaths())); source('$d')" > $S/c3_$(basename $d).log 2>&1 && echo "$d OK $(grep -c PASS $S/c3_$(basename $d).log) PASS" || echo "$d FAIL"; done > bench/proto/r12/memory/tests_c3.log
bench/proto/r12/memory/run_c3.sh > bench/proto/r12/memory/run_c3.log 2>&1
