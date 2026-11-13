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

## Round 5 — step 1 Stan model (memory specialist's mem_stan), evaluated separately
Same model and objective (log density at the baseline optimum identical to 1e-10); fewer autodiff nodes
(one matrix-vector product for the fixed components, one tcrossprod for the factors). NOT bit-identical:
L-BFGS takes a different path. bench/stan_downstream.R (learn_covariance P = 10, then step 2 on 40 genes
with each kernel):
| data | step-1 time base / new | cov_suez max rel diff | max d(-log10 p_geno) | max d(-log10 p_interact) | lead SNP same | Spearman |
| S | 1.61 / 1.06 s | 1.4e-4 | 0.0018 | 0.0013 | 100% | 0.9999996 |
| M | 7.53 / 4.27 s | 3.1e-4 | 0.0025 | 0.0026 | 100% | 0.9999992 |
Peak RAM (memory specialist): S 372 -> 326 MB, M 686 -> 431 MB, L 2.77 -> 1.12 GB.
tests/learn_covariance_checks.R: 6/6 PASS. Accepted as a numerical (not methodological) change, documented.

## Round 6 — combined candidate p9 (numerical specialist, reviewer's combination), integrated
p9_combined.cpp -> src/lmm.cpp (suez_fit_gene and suez_fit_models on one engine). Kernel level (interleaved):
1.9x (M_eqtl) .. 8.1x (M_fac3) vs p0, ~0.83x of p7's speed; mem_covs 325 MB vs baseline 326 MB (p7 2.6 GB);
extreme_cases <= 1.8e-10 from exact (baseline 2.7e-10); LA fit_models ll 1.7e-13, 7.4x.
End-to-end (bench/validate_all.sh, bench/results/validate_r6_p9.txt): all 18 datasets EQUIVALENT
(data ll <= 1e-9, boot ll <= 1e-5, beta/se 1e-6, identical NA/df, identical lead SNPs and significant sets);
tests/ 30/30 PASS. Scouting wall (1 core, 1 rep): S 4.8x, M 5.6x, L 4.9x, XL 2.15x, M_grm 3.0x, M_learn 4.9x,
M_num 2.8x, M_fac3 5.6x, M_eqtl 1.6x, M_cov 6.5x, M_miss 2.6x, M_dense 6.9x, S_rare 6.7x, LA 5.4x;
peak RSS 0.59-0.94x. Authoritative repeated benchmarks: final round.

## Round 7 — kernel_contrast_eigen at large N (bottleneck after round 6: 95% of XL, 60% of L)
Profile (bench/profile_e2e.R, r6_p9): XL eigen() 12.7 s + O(n^3) Q'KQ and QV products 7.4 s.
Prototypes (bench/proto/round7/eig_bench.R, GRM-like K): base 2.58 s (n = 1000) / 23.1 s (n = 2000);
Householder projection 0.01/0.02 s; R eigen() on HKH 0.65/13.2 s; Eigen SelfAdjointEigenSolver 0.26/2.48 s.
Eigenvalues agree to 7e-13; projectors on distinct eigenvalues to 2e-15.
Integrated (Householder + Eigen solver, .Call suez_eigen_sym). Consequence: inside tied eigenvalues the
basis differs, so the bootstrap draws for a given seed change (data fits are basis-invariant).
Validation (bench/results/validate_r7_eig.txt, validate_boot_r7_eig.txt): data l <= 8e-12, beta/se,
NA/df, lead SNPs, significant sets and eGenes identical on all 18 datasets; bootstrap null LR distribution
vs baseline KS p 0.11-0.98 (M_unbal 0.012; 17 tests), quantiles to q99.9 matching. tests/ 30/30 PASS.
Wall vs r6 (1 rep): XL 2.68x, L 1.27x, L_grm 1.22x, M_fac3 1.27x, M_unbal 1.25x, others 1.0-1.14x.
Accepted as a documented reproducibility change (same bootstrap distribution, different draws per seed).
History cleanup: sourceCpp .cache artifacts (475 MB) removed from the branch's commits.

## Round 8 — reviewer's required fixes to p9 (final review, bench/proto/review/cand/)
1. COLLINEAR pivot tolerance 1e-12 -> 1e-14 (collinear_sweep.R: 1e-12 dropped real information at
   1 - R^2 = 1e-12, l off by 0.019; 1e-14 exact there; still drops at 1 - R^2 = 1e-14, declared as a numerical
   rank tolerance for degenerate designs).
2. beta/se on the columns the pivoted LDLT keeps, NA for dropped columns (0.4 and p9 inverted a singular
   X'WX: NaN, or se = 0 / p = 0 in p9 for identified local-ancestry terms when an ancestry is absent).
   rankdef_beta.R: beta/se = identified fit to 1e-11 (no NaN); la_edge.R case D: main effects equal 0.4, the
   interaction SEs 0.4 reported as NaN are now finite; all LA cases data l <= 9.4e-13.
3. Comments say "suez 0.4" instead of "p0".
New regression tests: kernel x1e6 / x1e-6 invariance (model_checks.R), absent third ancestry equals the
two-ancestry fit (local_ancestry_checks.R; fails on 0.4, passes now). tests/: 32 PASS.
r8 vs r7: EQUIVALENT on all 18 datasets incl. bootstrap (validate_r8_vs_r7.txt); vs baseline: data columns
within tolerance, 0 NA mismatches (validate_r8_fix.txt).
