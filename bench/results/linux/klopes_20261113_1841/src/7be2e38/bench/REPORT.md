# suez performance optimization: final report

Branch `perf-optimization` (suez 0.5) vs git tag `perf-baseline` (suez 0.4). Apple M5 (4 performance + 6
efficiency cores), 16 GB, macOS, R 4.6.0 with the reference BLAS/LAPACK, clang -O2. Full chronology:
`bench/EXPERIMENT_LOG.md`; specialist reports: `bench/proto/<role>/REPORT.md`.

## 1. Headline

Median of 3 serial runs, baseline and final interleaved in each repetition; every individual repetition
improves (`speedup_min` > 1 everywhere). Raw rows: `bench/results/runs_final.tsv`; table:
`bench/results/final_summary.tsv`; figures: `bench/results/final_*.png`.

| workload | N | tests | 0.4 wall (s) | 0.5 wall (s) | speed-up [min-max of reps] | peak RSS 0.4 -> 0.5 (MB) |
|---|---|---|---|---|---|---|
| S (default kernel) | 200 | 19,942 | 3.48 | 0.62 | 5.6x [5.5-5.7] | 322 -> 193 (-40%) |
| M | 400 | 60,199 | 18.85 | 2.80 | 6.7x [6.4-6.8] | 362 -> 224 (-38%) |
| L | 1000 | 59,949 | 49.42 | 7.48 | 6.6x [6.6-6.6] | 527 -> 412 (-22%) |
| XL | 2000 | 23,859 | 62.94 | 10.63 | 5.9x [5.9-6.0] | 984 -> 736 (-25%) |
| M_grm (GRM kernel) | 400 | 30,049 | 9.40 | 2.93 | 3.2x | 345 -> 215 (-38%) |
| L_grm (GRM kernel) | 1000 | 39,899 | 33.40 | 10.83 | 3.1x | 512 -> 368 (-28%) |
| M_learn (learned kernel) | 400 | 30,049 | 8.76 | 1.45 | 6.0x | 344 -> 206 (-40%) |
| M_lc (learn_covariance kernel) | 400 | 30,049 | 9.31 | 1.55 | 6.0x | 346 -> 206 (-40%) |
| M_num (numeric condition, GRM) | 400 | 30,049 | 9.40 | 2.97 | 3.2x | 345 -> 215 (-38%) |
| M_fac3 (3 levels) | 600 | 30,049 | 21.31 | 2.94 | 7.2x | 372 -> 243 (-35%) |
| M_eqtl (no interaction, GRM) | 200 | 30,049 | 1.88 | 0.98 | 1.9x | 326 -> 193 (-41%) |
| M_cov (10 covariates) | 400 | 30,049 | 27.56 | 3.75 | 7.4x | 338 -> 218 (-36%) |
| M_miss (2% missing genotypes) | 400 | 30,049 | 12.42 | 2.73 | 4.5x | 345 -> 207 (-40%) |
| M_dense (overlapping windows) | 400 | 148,125 | 46.09 | 6.50 | 7.1x | 517 -> 388 (-25%) |
| M_unbal (1-4 samples/individual) | 488 | 30,049 | 11.52 | 1.67 | 6.9x | 346 -> 210 (-39%) |
| S_rare (rare SNPs, strong LD) | 200 | 19,104 | 3.27 | 0.59 | 5.5x | 322 -> 193 (-40%) |
| HS (high signal, GRM) | 400 | 19,950 | 6.29 | 2.02 | 3.1x | 326 -> 210 (-36%) |
| LA (local ancestry) | 300 | 9,950 | 9.76 | 1.88 | 5.2x | 367 -> 225 (-39%) |
| step 1, S (learn_covariance P = 10) | 200 | - | 1.77 | 1.13 | 1.6x | 380 -> 339 (-11%) |
| step 1, M | 400 | - | 8.02 | 4.12 | 1.9x | 692 -> 470 (-32%) |

Multi-core (doMC forks; peak RSS summed over the process tree, an upper bound):

| dataset | cores | 0.4 wall (s) | 0.5 wall (s) | speed-up | tree RSS 0.4 -> 0.5 (MB) | parallel efficiency 0.4 / 0.5 |
|---|---|---|---|---|---|---|
| M | 2 | 10.52 | 1.60 | 6.6x | 822 -> 522 | 0.90 / 0.87 |
| M | 4 | 6.40 | 1.00 | 6.4x | 1240 -> 888 | 0.74 / 0.70 |
| M | 8 | 4.49 | 0.69 | 6.5x | 2157 -> 1393 | 0.52 / 0.51 |
| M_dense | 2 | 27.68 | 3.80 | 7.3x | 983 -> 676 | 0.83 / 0.86 |
| M_dense | 4 | 18.94 | 2.37 | 8.0x | 1484 -> 1076 | 0.61 / 0.68 |
| M_dense | 8 | 14.39 | 1.62 | 8.9x | 2358 -> 1726 | 0.40 / 0.50 |

Disk I/O: 0 block reads/writes in every run (inputs are in-memory RDS; checkpoint I/O measured separately by
the memory specialist: no measurable write cost, resume 40 ms/file).

## 2. Baseline and bottlenecks

Round-0 profile of 0.4 (`bench/proto/p0_timed.cpp`): the grid search over h was 60-65% of C++ time (one
data grid ~20%, three bootstrap grids ~40%; each grid point allocated n x S temporaries), Newton refinement
~30%, rotation U'G 2-11%. Loading the package cost 220 MB RSS (rstan attached at load), per forked worker too.
After the grid was fixed, the bottleneck moved to `kernel_contrast_eigen` at large N (95% of XL: `eigen()`
with the reference LAPACK and O(N^3) dense products), then to rotation (17-71% of C++ time) and, with missing
genotypes, to `easy_impute`'s `svd()` (51% of M_miss).

## 3. Specialist agents

| specialist | main findings | outcome |
|---|---|---|
| Algorithm | eigenvalue-group compression with per-group QR (p3); Householder contrast basis; rotation caching across overlapping windows | compression idea adopted (via p9); Householder adopted in round 7; rotation caching prototyped in round 9 and shelved |
| Numerical | grid as one GEMM, Schur solves, group moments (p7, fastest), duplicate-row rotation, no copies; -O3/-mcpu=native and log-det tricks rejected | p7 engine + reviewer fixes = p9, integrated |
| Memory | rstan attached at load = 150-185 MB; N x N temporaries; quantile_normalize 6 copies; doMC GC breaks copy-on-write; step 1 autodiff ~ (2P + 3P_adj + 5) N^2 | lazy rstan, lean set-up, Stan rewrite adopted |
| Parallel | efficiency loss at 8 cores = slower E cores (2.4x), not imbalance; dynamic scheduling costs 14-33 ms/gene; threads 3.1x in-process but slower than forks end-to-end | kept static doMC; threads not adopted (opt-in memory option documented) |
| Reviewer | harness audit (11 issues incl. multi-core RAM, biased kernels); validation and edge-case suites; found bugs in p3 (-Inf grid, l off by 175), p7 (null-eigenvalue ties, l off by 0.026), p9 (collinearity tolerance, singular Wald SEs) and in round 10 (genes silently dropped) | every bug fixed before integration; final verdict: stable and maintainable |

## 4. Prototypes and leaderboard

Kernel-level speed-up vs 0.4's C++ (interleaved, median of 3; `bench/results/h2h_round1.tsv` and agent
reports), correctness vs 0.4.

| rank | prototype | strategy | speed-up (range) | RAM | correctness | complexity | status |
|---|---|---|---|---|---|---|---|
| 1 | p9_combined (numerical + reviewer) | eigen-group moments, null-group Taylor, SNP blocks | 1.9-8.1x | = 0.4 with 30 covariates | data l <= 1.4e-12; extreme cases <= 2e-10 from exact | high (670 lines, documented) | **integrated** (+ round 8 fixes) |
| 2 | p7_exact (numerical) | group moments, Schur, dedup rotation | 2.3-9.3x | 2.1-2.6 GB with 30 covs | l off by 0.026 near h -> 1 | medium | superseded by p9 |
| 3 | p3_compress_fast (algorithm) | per-group QR compression | 1.8-5.7x | 2.4 GB with 30 covs | l off by 175 at kernel x1e6 | high | superseded |
| 4 | p2_gemm (numerical) | grid as GEMM + Schur | 2.2-2.8x | - | exact | medium | superseded |
| 5 | p7_tol0 (reviewer) | p7 with exact ties only | 3.9-7.3x | - | exact | medium | superseded |
| 6 | p4_dedup / p5_loops | rotation over repeated rows; plain-loop solves | +20% rotation; +6-15% | - | exact | low | folded into p9 |
| 7 | p1_profile | no copies / allocations | 1.1x | - | exact | low | folded into p9 |
| - | p6_profile, p8_logdet, -O3, -mcpu=native | combined Profile product; frexp log-det; compiler flags | 0.95-1.09x, 0.67x, <= 1.03x | - | exact | - | rejected |
| - | Householder + Eigen solver (round 7) | O(N^2) basis, Eigen eigen-solver | 9x on the set-up at N = 2000 | - | data exact; new bootstrap draws | low | **integrated** |
| - | easy_impute via Gram eigen (round 10) | rank-k projection | 2.2x imputation | - | identical genotypes | low | **integrated** |
| - | rotation caching (round 9) | chunks of genes, fixed 64-SNP rotation blocks | 1.16-1.20x at 1 core; 0.91x at 4 cores (M) | +16..155 MB | identical across cores | high (+120 lines) | **rejected / shelved** |
| - | mem_lazy, mem_lean, mem_lean2 | lazy rstan, lean set-up | - | -30..-40% | identical() | low | **integrated** |
| - | mem_stan | Stan model with fewer autodiff nodes | 1.6-2.0x step 1 | -11..-60% step 1 | same objective, ~1e-4 path difference | low | **integrated** (documented) |
| - | dynamic doMC, chunk4/8 | scheduling | 0.79-1.00x, 0.84-1.12x | - | identical | low | rejected |
| - | std::thread over grid/SNPs | threads | 3.1x in-process; slower than forks | 0.35 vs 2.1 GB at 8 workers | identical | medium | not adopted (option) |

## 5. Statistical validation

- **Tolerances** (fixed after round 1, `bench/proto/review/validate.R`): data log-likelihoods 1e-9; beta and
  se relative 1e-6 (beta relative to max(|beta|, se)); identical NA pattern and df. The bootstrap columns
  were at first required to match within 1e-5, since 0.4 itself moves them by 7.6e-7 under a 1e-15
  perturbation of y. After round 7 they are judged by their distribution instead (see below).
- **End-to-end, 18 datasets** (`bench/results/equivalence_final2.tsv`): max |dl| 4.4e-11, max
  |d -log10 p| 9.1e-8, max |dbeta|/se 7.3e-8, max relative dse 2.7e-8, 0 NA mismatches, identical lead SNP
  per gene, identical significant sets (p < 1e-3, 1e-5, 1e-8) and eGene sets (Bonferroni + BH 5%), identical
  power/FDR vs the simulated truth.
- **Bootstrap**: same distribution as 0.4 (`bench/results/validate_boot_final2.txt`): KS p 0.11-0.98 on 16
  datasets, 0.012 on M_unbal (17 tests); quantiles to q99.9 match. On a null run with permuted genotypes
  (reviewer), lambda vs chi2 is 1.057 (0.4) vs 1.060 (0.5) on S and 0.997 vs 1.037 on M_grm.
- **Null calibration** (permuted genotypes, same seed): data columns identical to 0.4; lambda_GC
  1.03/1.01 (S p_geno/p_interact) and 0.98/0.98 (M_grm), same as 0.4.
- **Exactness**: extreme cases (h -> 1 with null directions, kernels x1e6/x1e-6, spread null eigenvalues)
  within 2.2e-10 of an exact optimiser (0.4: 2.7e-10); bimodal and random-kernel stress tests: 0 misses;
  near-perfect fits closer to exact than 0.4 (1.5e-7 vs 0.023).
- **Reproducibility**: identical() across 3 repetitions and across 1/2/4/8 cores (M, M_dense).
- **Tests**: `tests/` 33 checks PASS (3 new regression tests); `R CMD check --no-manual` (vignette rebuilt):
  Status 1 WARNING, the environment's "Rcpp built under R 4.6.1", as for 0.4.

**Changes that are not bit-for-bit equivalent** (documented in CHANGELOG 0.5, none changes the statistical
model):
1. Bootstrap draws for a given seed differ from 0.4 (eigenvectors inside tied eigenvalues are an arbitrary
   basis; each noise component is independent N(0, s2k lambda_i + s2), so the distribution is invariant).
2. Step 1 converges along a different optimiser path: `cov_suez` ~1e-4 relative; step-2 -log10 p within 0.003,
   identical lead SNPs (`bench/stan_downstream.R`).
3. Degenerate designs: a column with 1 - R^2 < 1e-14 on the others is treated as absent and its Wald estimate
   is NA (0.4 inverted the singular matrix and gave NaN or arbitrary SEs).

## 6. Final implementation changes

- `src/lmm.cpp`: new fitting engine (eigen-group moments, GEMM grid, Newton on moments with direct-residual
  fallback, SNP blocks, duplicate-row rotation, Wald on identified columns), shared by `suez_fit_gene` and
  `suez_fit_models`; `suez_eigen_sym`.
- `R/lmm.R`: Householder contrast basis + Eigen solver; column-wise symmetry check.
- `R/map_interaction_qtl.R`: Z Z' default kernel freed after use; row-wise `quantile_normalize`; Gram-eigen
  `easy_impute`; exact non-constant filter.
- `src/stan_files/suez_step_1.stan`: fewer autodiff nodes.
- `NAMESPACE` / `R/suez-package.R` / `R/learn_covariance.R`: rstan not attached at load.
- `tests/`: kernel-scale invariance, absent third ancestry, monomorphic SNP with missing calls.
- Public interface unchanged; no new dependencies.

## 7. Stopping assessment

- Optimization rounds: 10 (grid/refine; parallel; memory study; memory integration; step 1; combined engine;
  eigen set-up; reviewer fixes; rotation caching; R-side costs).
- Distinct prototypes: 26 (listed above and in the agent reports), 4 strategy families (algorithmic
  reformulation, numerical implementation, memory, parallel execution).
- Best runtime improvement: 7.4x on one core (M_cov; 86% less wall time), 8.9x at 8 cores (M_dense);
  representative S/M/L/XL: 5.6-6.7x.
- Best RAM improvement: 41% less peak RSS (M_eqtl); 22-41% on every step-2 workload; 32% less in step 1 (M);
  27-35% less process-tree RAM at 2-8 cores.
- Statistical validation: passed (all 18 datasets, edge/extreme/null suites, tests, R CMD check).
- Diminishing returns: round 9 (the last high-potential idea, rotation caching) gave only 14-17% at one core
  for large N, a 10% regression at 4 cores and more RAM, so it was rejected. Round 10 only gained much
  where imputation dominates (M_miss 1.65x; M_dense 1.08x; other workloads unchanged). The 10-round
  budget is reached.
  The last two rounds are not both below 5% on every workload, so criterion E.1 is not strictly met; the stop
  rests on F: the remaining identified opportunity needs disproportionate complexity for a mixed result.
- Remaining opportunities (not pursued): rotation caching as an opt-in for large-N single-core runs (+14-17%);
  opt-in C++ threads for memory-bound machines (1/6 of the RAM of 8 forks at ~75% of their speed); a streaming
  result writer for genome-wide runs (the result table dominates RAM above ~10M rows; an interface change);
  an optimised BLAS/LAPACK in the user's R (outside the package).

## 8. Remaining limitations

- GRM-like kernels with all-distinct eigenvalues gain 3.1-3.2x (1.9x without interaction terms) rather than
  6-7x: the grouping has nothing to compress, and rotation plus Newton on n contrasts remain.
- Rotation U'G is now the largest single cost at large N (53-71% of fitting time on L/XL).
- Benchmarks are simulated; real data with highly uneven gene sizes were not available.
- Timings are on one Apple M5; the reference BLAS was used throughout (as in a default R install).

## 9. Reproduce

```bash
git checkout perf-optimization
mkdir -p bench/lib/baseline bench/lib/final
git worktree add /tmp/suez-baseline perf-baseline && R CMD INSTALL --no-docs --library=bench/lib/baseline /tmp/suez-baseline
R CMD INSTALL --no-docs --library=bench/lib/final .
Rscript bench/sim.R                         # datasets -> bench/data/ (~400 MB)
bench/final_bench.sh 3                      # serial benchmark, ~40 min -> bench/results/runs_final.tsv
Rscript bench/analyze_final.R               # tables and figures in bench/results/
bench/validate_boot.sh final                # end-to-end validation vs baseline (data + bootstrap distribution)
Rscript bench/summarize_equivalence.R final # -> bench/results/equivalence_final.tsv
Rscript bench/proto/review/edge_cases.R lib:bench/lib/final
for t in tests/*.R; do Rscript -e ".libPaths(c('bench/lib/final', .libPaths())); source('$t')"; done
R CMD build . && R CMD check --no-manual suez_0.5.tar.gz
```
Kernel-level comparisons: `Rscript bench/proto/h2h.R 20 S M_grm L` (prototypes listed in its `CANDS`).
