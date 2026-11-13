# Linux checks

suez 0.5 was optimised and benchmarked only on macOS (Apple clang, R's reference BLAS/LAPACK). The Linux
machines that run it (GCC, libstdc++, OpenBLAS, no GPU) differ in exactly the places round 11 found problems
(`bench/EXPERIMENT_LOG.md`, round 11). `run_checks.sh` checks them on the machine it runs on.

## Run

From the repository root, with the R that will use suez on `PATH` (it needs rstan, StanHeaders, RcppEigen,
RcppParallel, BH, doMC, dplyr, readr, magrittr, foreach and a C++17 compiler):

```bash
bash bench/linux/run_checks.sh
```

With the eqtl_smk pixi environment:

```bash
pixi run --manifest-path /path/to/PARDoS_GWAS/eqtl_smk/tests/pixi.toml bash bench/linux/run_checks.sh
```

Tests only (no timings, a few minutes):

```bash
bash bench/linux/run_checks.sh --quick
```

Options (environment variables): `REFS="d012046 7be2e38"` (refs timed against the working tree: 0.3 and
0.5), `NS="1000 2000 5000"`, `REPS=3`, `OUT=<dir>`, `JOBS=<make jobs>`. The full run takes about 30-60 min
(builds ~2 min each; the N = 5000 timings dominate). Keep the machine otherwise idle for the timings.

To also time against R's reference LAPACK on a FlexiBLAS system, run it a second time with
`FLEXIBLAS=netlib OUT=bench/results/linux/<host>_netlib`.

## What it does

1. Records the environment (`environment.txt`: compiler, flags, R, BLAS/LAPACK with symlinks resolved).
2. Builds the working tree (`cur`), and a second copy with `-O1 -Wall -Wextra -D_GLIBCXX_ASSERTIONS`
   (`assert`: libstdc++ bounds and comparator checks, so undefined behaviour in `std::sort` or vector access
   aborts instead of passing silently). Compiler warnings from `src/lmm.cpp` go to `lmm_warnings.txt`.
3. Runs `tests/*.R` against both builds (the log shows which eigen solver `auto` chose).
4. Times `kernel_contrast_eigen`'s solvers (`bench/proto/round11/eig_bench.R`) with 1 and the default number
   of BLAS threads: `eig_threads_*.tsv` (eigen_cpp = Eigen's, lapack = `eigen()`).
5. Builds `REFS` and times `suez_fit_gene` on identical inputs (20 genes x 50 SNPs, GRM-like kernel, with and
   without an interaction, 1 BLAS thread, interleaved): `fit_summary.tsv` has the median time, the speed-up
   vs the first ref and the largest difference of every output column vs that ref.

All output goes to `bench/results/linux/<host>_<date>/` (`summary.txt` first).

## What to look for

- Any `FAIL`, `ERROR` or a failed build, especially in the `assert` build only (undefined behaviour that the
  normal build hides).
- `eig_threads_*.tsv`: if `eigen_cpp` beats `lapack` where `auto` chose `lapack` (or the reverse), the
  detection in `R/lmm.R`'s `eigen_solver()` is wrong for this BLAS.
- `fit_summary.tsv`: `cur` should not be slower than 0.3 (`d012046`) at any N; `max_abs_diff_vs_ref` should
  stay below ~1e-8 with `na_mismatch` 0.
