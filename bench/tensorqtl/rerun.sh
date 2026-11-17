#!/usr/bin/env bash
# New benchmark round with the current suez (machine otherwise idle; pull first, so this script is current):
#   git pull && bash bench/tensorqtl/rerun.sh         # re-time suez only, keep the tensorQTL results (~30 min)
#   git pull && bash bench/tensorqtl/rerun.sh all     # re-time everything, tensorQTL included (~1 h)
# Steps: build suez at HEAD into bench/lib/cur; run tests/ against that build (stops on a failure);
# simulate / export any missing dataset; move the previous outputs to bench/tensorqtl/out_<suez>_<date>/
# (git-ignored; suez mode keeps the tensorQTL runs in out/); run_all.sh; analyze.py (results/ and
# doc/suez_vs_tensorqtl.html). Logs of the build and tests: bench/lib/.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
mode=${1:-suez}
case $mode in suez|all) ;; *) echo "usage: bash bench/tensorqtl/rerun.sh [suez|all]"; exit 1 ;; esac
P() { pixi run --manifest-path bench/tensorqtl/pixi.toml "$@"; }
T=bench/tensorqtl; O=$T/out
step() { echo; echo "== $(date +%H:%M:%S) $*"; }

rev=$(git rev-parse --short HEAD)
[ -z "$(git status --porcelain -- R src DESCRIPTION NAMESPACE)" ] || echo "note: local changes in R/ or src/ are built too"

step "build suez $rev into bench/lib/cur"
mkdir -p bench/lib/cur
P -e r R CMD INSTALL --no-docs --library=bench/lib/cur . > bench/lib/install.log 2>&1 ||
  { tail -30 bench/lib/install.log; echo "build failed (bench/lib/install.log)"; exit 1; }
tail -1 bench/lib/install.log

step "package tests"
for t in model_checks local_ancestry_checks learn_covariance_checks; do
  R_LIBS=bench/lib/cur P -e r Rscript tests/$t.R > bench/lib/test_$t.log 2>&1 ||
    { grep -E '^FAIL|Error' bench/lib/test_$t.log; echo "tests failed (bench/lib/test_$t.log)"; exit 1; }
  echo "$t: $(grep -c '^PASS' bench/lib/test_$t.log) checks passed"
done

step "datasets"
DS="S M L XL XXL G M_grm L_grm XL_grm M_fac3 M_num M_cov M_eqtl M_unbal HS NULLR REAL"
for ds in $DS; do [ -f bench/data/$ds.rds ] || { P -e r Rscript $T/sim_tq.R; break; }; done
todo=""; for ds in $DS; do [ -f $T/data/$ds/geno.i8 ] || todo="$todo $ds"; done
[ -z "$todo" ] || P -e r Rscript $T/export.R $todo
echo "datasets ready"

step "previous outputs"
if [ -d $O ]; then
  old=$(sed -n 's/^suez: \([0-9a-f]*\).*/\1/p' $T/results/environment.txt 2>/dev/null); old=${old:-prev}
  B=${T}/out_${old}_$(date +%Y%m%d_%H%M)
  if [ $mode = all ]; then
    mv $O $B
  else
    cp -r $O $B
    rm -f $O/speed/*_suez* $O/stats/*_suez*
    [ -f $O/runs.tsv ] && awk -F'\t' 'NR == 1 || $3 !~ /^suez/' $O/runs.tsv > $O/runs.tmp && mv $O/runs.tmp $O/runs.tsv
  fi
  echo "previous outputs (suez $old) copied to $B (git-ignored)"
else
  echo "none"
fi

step "benchmark ($mode)"
bash $T/run_all.sh

step "report"
P -e py python $T/analyze.py

step "done: suez $rev"
column -t $T/results/speed_wide.tsv
git status --short $T/results doc/suez_vs_tensorqtl.html
echo
echo "To publish: git add $T/results doc/suez_vs_tensorqtl.html && git commit -m 'bench results on linux ($rev)' && git push"
