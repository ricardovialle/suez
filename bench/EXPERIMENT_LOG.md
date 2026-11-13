# Experiment log

Hardware: Apple M5 (4 performance + 6 efficiency cores), 16 GB RAM, macOS; R 4.6.0 with the reference
BLAS/LAPACK; clang++ -O2 (R defaults). Baseline: git tag `perf-baseline` (suez 0.4), installed in
`bench/lib/baseline`.

## Round 0 — baseline and profile

- Harness: `bench/sim.R` + `bench/datasets.R` (14 workloads), `bench/run.R` + `bench/bench.sh`
  (/usr/bin/time -l: wall, CPU, max RSS, I/O blocks), `bench/proto/harness.R` (per-gene kernel inputs).
- C++ time split (bench/proto/p0_timed.cpp): grid evaluation 60-65% (data grid ~20%, three bootstrap
  grids ~40%), Newton refinement ~30%, rotation U'G 2-11% (grows with n).
- Fixed RAM: library(suez) = 290 MB RSS vs 70 MB for bare R (rstan/StanHeaders loaded at load time).
- Specialist agents launched: algorithm, numerical, memory, parallel, independent reviewer
  (briefs: bench/CONTEXT.md; reports: bench/proto/<role>/REPORT.md).

## Round 0b — harness fixes after reviewer audit
- bench.sh: CPU from /usr/bin/time (forks included), process-tree peak RSS (memwatch.sh), thread caps.
- New datasets: L_grm, M_unbal, M_lc (kernel from learn_covariance), HS (high signal).
- Scouting baseline (1 core, wall): S 3.3 s, M 18.2, L 49.3, XL 62.4, M_grm 9.4, M_learn 8.8, M_num 9.5,
  M_fac3 21.5, M_eqtl 1.9, M_cov 28.2, M_miss 12.5, M_dense 52.9, S_rare 4.8, LA 11.7; step1 S 2.2 s.

## Round 1 — grid/refine bottleneck (algorithm + numerical specialists)
Prototypes: alg p2_compress, p3_compress_fast (eigen-group QR compression, direct residuals);
num p1_profile (no copies/allocs, 1.1x), p2_gemm (grid as GEMM + Schur solves, 2.2-2.8x),
p3_grouped (eigen-group moments), p4_dedup (rotation over repeated sample rows), p5_loops, p6_profile
(rejected, -5..+9%), p7_exact (p5 + residual fallback + LU se), p8_logdet (rejected, 0.67x);
-O3 / -mcpu=native rejected (<= 3%). Householder contrast eigen: saves 21% on XL but changes bootstrap
draws (basis inside repeated eigenvalues) -> held back.
Coordinator head-to-head (bench/proto/h2h.R, 20 genes, interleaved, CPU speed-up vs p0; bench/results/h2h_round1.tsv):
| dataset | p3 | p7 | max |dl| (both) | max |dl_boot| | dbeta/se |
| S | 5.70 | 7.82 | 1e-13 | 3.2e-7 | 6e-15 |
| M_grm | 2.87 | 3.98 | 2e-13 | 3.9e-7 | 2e-8 |
| M_lc | 5.39 | 8.20 | 8.5e-13 | 2.2e-7 | 2e-8 |
| M_eqtl | 1.83 | 2.29 | 7e-14 | 3.9e-7 | 4e-8 |
| M_cov | 3.75 | 9.33 | 2.3e-13 | 2.1e-7 | 9e-9 |
| L_grm | 2.83 | 3.83 | 9.7e-13 | 8.6e-7 | 4e-8 |
| L | 5.15 | 8.66 | 9.1e-13 | 5.6e-7 | 8e-15 |
No NA mismatches. Tolerance note: the baseline perturbed by 1e-15 relative in y moves l_boot by up to
7.6e-7 and beta by 4e-8 se (bench/proto/algorithm/floor.R), so the bootstrap columns are judged at 1e-6
with that floor in mind and beta by |dbeta|/max(|beta|, se) (reviewer H8).

## Round 2 — parallel execution (parallel specialist)
Prototypes: dynamic doMC (preschedule = FALSE), chunk4/chunk8 dynamic chunks, std::thread over grid points
and SNPs (opt-in suez.threads), 4 workers x 2 threads. Results identical() everywhere.
Verdict: dynamic rejected (fork cost 14-33 ms/gene, 0.79-1.00x); chunked no measurable gain (0.84-1.12x);
threads 3.1x in-process at 8 threads but slower than 8 forks end-to-end; their benefit is memory (0.35 GB
vs 2.1 GB on M at 8 workers). Efficiency loss at 8 cores is hardware (E cores 2.4x slower), not imbalance.
Key memory fact for round 3: each forked worker carries ~220 MB private, mostly rstan.

## Reviewer verdict on round 1 (bench/proto/review/cand/)
Both p3 and p7 "accept with changes"; no hidden methodological change (grid, peaks, Newton, ok_nested, RNG).
- p3 bug: grid marks valid points -Inf when exact-zero eigenvalues and the Taylor check fails at h < 1 (kernel
  scale >~1e5): l off by 175 (default x1e6), 119 (GRM). Fix p3_fix_anyzero.cpp (-4% speed).
- p7 bug: tie tolerance 1e-12*max groups null eigenvalues (1e-17..1e-12) at their mean -> model change near
  h -> 1 (l error 0.026 at x = 11.9; 0.003-0.008 with kernel x1e6). Tie-only grouping (p7_tol0) is exact
  (<= 5e-10 incl. weights 1e12) but loses ~half the speed on S. => normal-equation Profile + 1e-3 fallback is safe.
- Both: RAM blow-up with many covariates (30 covs, 1000 SNPs: p3 2.37 GB, p7 2.11 GB vs baseline 364 MB).
- Bootstrap tolerance: baseline vs itself (1e-15 perturbation) moves l_boot up to 7.6e-7; new acceptance rule
  |dl_boot| <= 1e-5, |d bootstrap LR| <= 1e-6, data l <= 1e-9 (tightened), identical NA/df.
- stress_exact (300 kernels) and new stress_bimodal: 0 misses for both.
Decision: build combined p9 (p7 engine + relative/null-aware grouping from p3 + cancellation guard + SNP
blocks + covariate pairs once) -> round 3.

## Round 3 (memory specialist) / Round 4 (integration of identical-output memory changes)
Prototypes: mem_lazy, mem_lean, mem_lean2, mem_stan (bench/proto/memory/REPORT.md).
Integrated mem_lean2 into the repo (lazy rstan, lean kernel_contrast_eigen/is_symmetric, quantile_normalize
row loop, Z Z' default kernel, U as Eigen::Map). tests/: 30 PASS, 0 FAIL. End-to-end output identical() to
baseline on S, M_grm, M_fac3, LA, M_lc. Peak RSS 1 core: S 323 -> 192 MB, M_lc 345 -> 203 MB.
mem_stan (step 1, -37% RAM, 1.9x faster on M, ~1e-4 relative change of the optimum) held for separate evaluation.
