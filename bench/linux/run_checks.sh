#!/usr/bin/env bash
# suez Linux checks (round 11): build, tests, a libstdc++-assertions build, eigen-solver and fit timings.
# Run from the repository root with the R that will use suez on PATH, e.g.
#   pixi run --manifest-path /path/to/eqtl_smk/tests/pixi.toml bash bench/linux/run_checks.sh
#   bash bench/linux/run_checks.sh --quick            # tests only (no timings), ~2 min after the builds
# Options (environment):
#   REFS="d012046 7be2e38"   git refs built and timed against the working tree (0.3 and 0.5 by default)
#   NS="1000 2000 5000"      sample sizes for the timings;  REPS=3  repetitions per timing
#   OUT=bench/results/linux/<host>_<date>   output directory;  JOBS=$(nproc)  make jobs per build
# Needs: R with the suez dependencies (rstan, StanHeaders, RcppEigen, RcppParallel, BH, doMC, dplyr, readr,
# magrittr, foreach), a C++17 compiler, git. Writes only under $OUT (builds go to $OUT/lib, sources to $OUT/src).
set -euo pipefail
QUICK=0; [ "${1:-}" = "--quick" ] && QUICK=1
ROOT=$(git rev-parse --show-toplevel); cd "$ROOT"
OUT=${OUT:-bench/results/linux/$(hostname -s)_$(date +%Y%m%d_%H%M)}
REFS=${REFS:-"d012046 7be2e38"}; NS=${NS:-"1000 2000 5000"}; REPS=${REPS:-3}
export MAKEFLAGS="-j${JOBS:-$(nproc 2>/dev/null || echo 4)}"
mkdir -p "$OUT/lib" "$OUT/src" "$OUT/logs" "$OUT/fit"
OUT=$(cd "$OUT" && pwd); R11="$ROOT/bench/proto/round11"
log() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$OUT/summary.txt"; }

# ---- environment ------------------------------------------------------------------------------------------
{
  echo "date: $(date)"; echo "host: $(uname -a)"; echo "suez HEAD: $(git rev-parse --short HEAD) dirty: $(git status --porcelain -- R src tests | wc -l | tr -d ' ') files"
  echo "cpu: $(lscpu 2>/dev/null | grep -m1 'Model name' | sed 's/.*: *//' || sysctl -n machdep.cpu.brand_string 2>/dev/null)"
  echo "cores: $(nproc 2>/dev/null || sysctl -n hw.ncpu)"; echo "CXX: $(R CMD config CXX17) / $($(R CMD config CXX17 | awk '{print $1}') --version | head -1)"
  echo "CXX17FLAGS: $(R CMD config CXX17FLAGS)"; echo "OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-unset} OMP_NUM_THREADS=${OMP_NUM_THREADS:-unset}"
  Rscript -e 'cat("R:", R.version.string, "\nLAPACK:", La_library(), "\nBLAS:", extSoftVersion()[["BLAS"]], "\nBLAS (resolved):", normalizePath(extSoftVersion()[["BLAS"]], mustWork = FALSE), "\n")' \
    -e 'for (p in c("rstan", "StanHeaders", "RcppEigen", "Rcpp", "RcppParallel", "BH")) cat(p, as.character(packageVersion(p)), "\n")'
} > "$OUT/environment.txt" 2>&1
cat "$OUT/environment.txt"

# ---- builds -----------------------------------------------------------------------------------------------
build() {   # build <label> <source dir> [makevars file]
  local lab=$1 src=$2 mv=${3:-}
  mkdir -p "$OUT/lib/$lab"
  if R_MAKEVARS_USER=${mv:-/dev/null} R CMD INSTALL --preclean --no-docs --library="$OUT/lib/$lab" "$src" > "$OUT/logs/install_$lab.log" 2>&1; then
    log "build $lab: ok ($(grep -c 'warning:' "$OUT/logs/install_$lab.log" || true) compiler warnings, $(grep 'warning:' "$OUT/logs/install_$lab.log" | grep -c 'lmm.cpp' || true) in lmm.cpp)"
  else log "build $lab: FAILED (see logs/install_$lab.log)"; tail -20 "$OUT/logs/install_$lab.log"; return 1; fi
}
rm -rf "$OUT/src/cur"; mkdir -p "$OUT/src/cur"
git ls-files -z --cached --others --exclude-standard -- . ':!bench' | tar --null -T - -cf - | tar -xf - -C "$OUT/src/cur"
build cur "$OUT/src/cur"
# libstdc++ debug assertions (bounds checks, std::sort comparator checks) + warnings, for the tests only
printf 'CXX17FLAGS = -O1 -g -Wall -Wextra -Wno-unused-parameter -Wno-ignored-attributes -D_GLIBCXX_ASSERTIONS\n' > "$OUT/Makevars.assert"
rm -rf "$OUT/src/assert"; cp -R "$OUT/src/cur" "$OUT/src/assert"
build assert "$OUT/src/assert" "$OUT/Makevars.assert" || true
grep 'lmm.cpp.*warning:' "$OUT/logs/install_assert.log" | sort -u > "$OUT/lmm_warnings.txt" || true

# ---- tests (each build) -----------------------------------------------------------------------------------
for lab in cur assert; do
  [ -d "$OUT/lib/$lab/suez" ] || continue
  for t in tests/*.R; do
    n=$(basename "$t" .R)
    if R_LIBS="$OUT/lib/$lab" Rscript -e "stopifnot(startsWith(find.package('suez'), '$OUT/lib/$lab'))" \
         -e "cat('eigen solver:', suez:::eigen_solver(), '\n')" -e "source('$t')" > "$OUT/logs/test_${lab}_$n.log" 2>&1; then st=ok; else st=ERROR; fi
    log "test $lab $n: $st, $(grep -c '^PASS' "$OUT/logs/test_${lab}_$n.log" || true) PASS, $(grep -c '^FAIL' "$OUT/logs/test_${lab}_$n.log" || true) FAIL ($(grep -m1 'eigen solver' "$OUT/logs/test_${lab}_$n.log" || true))"
  done
done
[ "$QUICK" = 1 ] && { log "quick mode: done"; exit 0; }

# ---- eigen solver timings (working tree build) ------------------------------------------------------------
for th in 1 default; do
  if [ "$th" = 1 ]; then export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1; else unset OPENBLAS_NUM_THREADS OMP_NUM_THREADS; fi
  Rscript "$R11/eig_bench.R" "$OUT/lib/cur" "$OUT/eig_threads_$th.tsv" $NS > "$OUT/logs/eig_threads_$th.log" 2>&1 \
    && log "eigen timings (threads $th): $(cut -f1-4 "$OUT/eig_threads_$th.tsv" | tr '\n' ';')" || log "eigen timings (threads $th): FAILED"
done

# ---- suez_fit_gene timings: refs vs working tree, interleaved, 1 BLAS thread -------------------------------
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for ref in $REFS; do
  rm -rf "$OUT/src/$ref"; mkdir -p "$OUT/src/$ref"; git archive "$ref" | tar -x -C "$OUT/src/$ref"
  build "$ref" "$OUT/src/$ref" || true
done
for N in $NS; do for i in 0 1; do for r in $(seq "$REPS"); do for lab in $REFS cur; do
  [ -d "$OUT/lib/$lab/suez" ] || continue
  Rscript "$R11/fit_bench.R" "$OUT/lib/$lab" "$OUT/fit/${lab}_N${N}_i${i}_r${r}.rds" "$N" 20 50 "$i" >> "$OUT/logs/fit.log" 2>&1 \
    || log "fit $lab N=$N i=$i r=$r: FAILED"
done; done; done; done
first=$(echo "$REFS" | awk '{print $1}')
Rscript "$R11/compare_fit.R" "$OUT/fit" "$first" "$OUT/fit_summary.tsv" | tee -a "$OUT/summary.txt"
log "done: $OUT (summary.txt, environment.txt, eig_threads_*.tsv, fit_summary.tsv, logs/)"
