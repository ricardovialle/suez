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

## Round 9 (in progress, algorithm specialist) — rotation caching across overlapping cis windows
Rotation share of C++ time after round 8 (bench/proto/round8/rot_share.R): S 17%, M 28%, L 53%, L_grm 25%,
M_dense 29%, XL 71%.

## Round 10 — R-side per-gene costs (coordinator)
Profile (bench/profile_e2e.R, r8): easy_impute's svd() = 51% of M_miss; apply(G, 2, sd) = 10% of M_dense.
- easy_impute through the eigen-decomposition of the smaller Gram matrix (rank-k reconstruction = projection
  on the leading eigenvectors): 2.88 -> 1.30 s on all 150 M_miss windows, 0 of 119,501 imputed genotypes
  differ (bench/proto/round10/impute_bench.R); identical on wide/tall/tiny shapes.
- apply(G, 2, sd) > 0 -> exact non-constant test (colSums(G != first row)); same in local-ancestry mode.
- fix_diag (unused) removed.
End-to-end output identical() to r8 on M_miss, M_dense, S_rare, LA, M; tests pass.
Wall, median of 3 interleaved (bench/results/runs_round10.tsv): M_miss 4.48 -> 2.72 s (1.65x),
M_dense 6.95 -> 6.43 s (1.08x); peak RSS M_dense 365 -> 388 MB (recheck in the final benchmarks).

## Round 9 result — rotation caching (algorithm specialist): REJECTED / shelved
Implementation (bench/proto/algorithm/r9/r9_rot_on_f81c6c4.patch): rotation from individual-level genotypes
(W = rowsum(U, individual)), fixed 64-SNP rotation blocks (Eigen GEMM bits depend on matrix width, so fixed
blocks keep results independent of chunking/cores), chunks of consecutive genes as the foreach unit, per-gene
path kept for NA windows, permute, TabixFile and LA. Correct: data l <= 3.4e-13, bootstrap <= 7.5e-7 vs r8,
identical() across 1/2/4 cores, tests pass.
Timing vs an exact build of f81c6c4 (3 reps, median): L 0.845, XL 0.860, M_dense 0.835, M 0.985 at 1 core;
M 1.103 (slower) and M_dense 0.948 at 4 cores; RSS +16..+75 MB (1 core), +94..+155 MB (4-core tree).
Decision: gains 14-17% only for large-N / dense windows at 1 core, a regression at 4 cores on M, more RAM and
~120 lines of chunking logic -> disproportionate complexity; kept as a documented option.

## Final review and validation
Reviewer (final candidate): stable and maintainable; one blocking regression from round 10 — easy_impute
zeroed only missing entries, so a SNP constant among its calls (NaN after scale()) made the gene fail and be
silently dropped (debug = FALSE). Fixed (temp[is.na(temp)] = 0, as 0.4); regression test added (fails on the
buggy build, passes on 0.4 and 0.5); imputation identical to 0.4 on all of the reviewer's edge shapes.
Final serial benchmark (bench/final_bench.sh 3; bench/results/final_summary.tsv): step 2 1.9-7.4x at 1 core,
peak RSS -22..-41%; step 1 1.6-1.9x, -11..-32%; 8 cores 6.5x (M) / 8.9x (M_dense) with -27..-35% tree RAM.
The fix does not change any benchmark output (identical() on 9 datasets), so these timings stand.
End-to-end equivalence of the final build (bench/results/equivalence_final2.tsv): all 18 datasets PASS.
Reproducible: identical() across 3 runs and 1/2/4/8 cores. tests/ 33 PASS; R CMD check: environment WARNING only.
Report: bench/REPORT.md.

## Round 11 — review from a downstream user (eqtl_smk) on OpenBLAS; large N; Linux
Trigger: eqtl_smk's benchmark (suez 0.3 -> 0.5, its pixi R 4.5.3 with OpenBLAS) had identical results but its
5000-donor point slower (57 -> 65-103 s). Every earlier round used R's reference BLAS/LAPACK and clang on macOS;
the Linux machines run GCC/libstdc++ and OpenBLAS. Scripts and builds: bench/proto/round11/ (bench/lib/ob_*).
1. kernel_contrast_eigen's solver (round 7) depends on the LAPACK. GRM-like K, seconds
   (bench/results/round11_eig_*.tsv; lapack = same Householder projection + eigen()):
   | N | ref LAPACK: Eigen / lapack | OpenBLAS 1 thread: Eigen / lapack | OpenBLAS default threads: Eigen / lapack |
   | 1000 | 0.38 / 0.69 | 0.39 / 0.12 | 0.40 / 0.20 |
   | 2000 | 2.88 / 5.56 | 2.93 / 0.67 | 2.93 / 0.72 |
   | 5000 | - | 38.1 / 5.7 | 38.2 / 6.0 |
   Eigenvalues agree to 6e-13, leading projectors to 2e-13. Integrated: eigen_solver() (R/lmm.R) uses LAPACK when R's
   BLAS or LAPACK (symlinks resolved) is OpenBLAS/MKL/Accelerate/vecLib/BLIS/FlexiBLAS/ATLAS, else Eigen's;
   options(suez.eigen_solver) overrides. Unmeasured: reference LAPACK over a fast BLAS (e.g. macOS vecLib BLAS),
   classified as optimised. Bootstrap draws per seed change wherever LAPACK is now chosen (as in round 7).
2. suez_fit_gene slower than 0.3 at large N (bench/proto/round11/fit_bench.R, 20 genes x 50 SNPs, GRM-like K,
   OpenBLAS 1 thread, 3 interleaved reps; bench/results/round11_fit_summary_mac_openblas.tsv), speed-up vs 0.3:
   | N | cis 0.5 / r11 | interaction 0.5 / r11 |
   | 1000 | 1.3 / 1.4 (1 rep) | 1.4 / 1.5 (1 rep) |
   | 2000 | 1.00 / 1.21 | 1.07 / 1.44 |
   | 5000 | 0.89 / 1.37 | 0.79 / 1.45 |
   Cause: the genotypes were rotated inside the moment blocks, whose width is max(16, 4 MB / (8 R nper)) = 16 SNPs at
   N = 5000, so U (N x N-1, 200 MB) was streamed from memory once per 16 SNPs. Fix: rotate_all() rotates every
   column once per gene in fixed 256-column chunks (suez_fit_gene and suez_fit_models). Outputs vs 0.5 <= 3e-11,
   vs 0.3 <= 1.7e-9, identical NA patterns.
3. Linux-only risk: rotate() sorts sample rows with a comparator built on < and !=; a NaN in G or a term breaks
   strict weak ordering (undefined behaviour; libstdc++'s unguarded partition can read out of bounds, libc++ on
   macOS usually does not crash). map_interaction_qtl always imputes; direct .Call users (eqtl_smk) do not have to.
   Now refused with an error. New tests (model_checks.R): Eigen vs LAPACK solvers same log-likelihoods
   (3e-10 reference, 1.8e-9 OpenBLAS), NaN genotypes refused. tests/: 35 PASS with both R setups.
4. bench/linux/run_checks.sh (+ README): environment record, build, tests also under -D_GLIBCXX_ASSERTIONS
   -Wall -Wextra, solver timings at 1/default threads, fit timings vs 0.3/0.5. Dry-run on macOS OK; not yet run
   on Linux.
End-to-end, reference R (solver = Eigen, so this checks the rotation change): r11 vs 0.5 (final2) identical() on
17 of 18 datasets, L_grm bootstrap l within 6.3e-7; validate.R EQUIVALENT on all 18
(bench/results/validate_r11_vs_final2.txt). Vs the 0.4 baseline (validate_r11_ref.txt) the verdicts are the same
as 0.5's own (bootstrap draws, round 7).
Not checked here: GCC builds (no GCC on this Mac), R CMD check, end-to-end equivalence and the bootstrap
distribution with the LAPACK solver chosen (round 7 did both for the Eigen basis change).

# Round 12+ campaign (baseline = tag `perf-baseline-r12` = suez 0.5 + round 11, lib bench/lib/base12)
Brief: bench/CONTEXT_R12.md. Benchmark: bench/round12_bench.sh (interleaved implementations; RUNS tsv).
New workloads: XL_grm (N = 2000, GRM), XXL (N = 4000), M_cov30 (30 covariates), P_many (3000 genes x 10 SNPs).
Harness: bench/proto/harness.R with SUEZ_HARNESS_LIB=bench/lib/base12; bench/proto/r12/p12_base.cpp (base12
engine as sourceCpp, identical() to the installed build on 14 datasets), p12_timed.cpp + profile12.R (phases).

## Round 12.0 — baseline and profile (coordinator)
Scouting baseline, 1 rep (bench/results/runs_r12_baseline.tsv), 1 core wall / peak RSS: S 0.62 s / 194 MB,
M 2.92 / 231, L 7.79 / 423, XL 10.9 / 760, XXL 29.6 / 1264, M_grm 2.97, L_grm 10.7, XL_grm 12.1 / 656,
M_cov30 9.5, M_dense 6.5 / 402, P_many 6.8 / 247; XL_grm 2/4/8 cores 7.8 / 5.8 / 4.5 s (efficiency 0.33 at 8),
tree RSS 2.4 GB at 8. Step 1: S 1.11 s / 340 MB, M 4.06 s / 470 MB.
C++ phase shares (bench/results/r12_profile_base.tsv; % of fit_gene time):
| ds | N | rotate | data moments+grid | data refine+Wald | boot sim+moments | boot grid | boot refine |
| S | 200 | 14 | 16 | 16 | 15 | 14 | 19 |
| M | 400 | 29 | 14 | 11 | 20 | 9 | 14 |
| L | 1000 | 54 | 10 | 4 | 20 | 3 | 5 |
| XL | 2000 | 70 | 7 | 1 | 15 | 1 | 2 |
| XXL | 4000 | 86 | 3 | 0.4 | 9 | 0.3 | 0.5 |
| M_grm | 400 | 12 | 12 | 23 | 10 | 11 | 30 |
| XL_grm | 2000 | 37 | 8 | 16 | 8 | 7 | 23 |
| M_cov30 | 400 | 2 | 25 | 18 | 8 | 21 | 23 |
R level (bench/results/r12_profile_e2e_base.log): XXL kernel_contrast_eigen 74% of wall (16.7 s); XL_grm 19%;
P_many: cat() 38% of wall (one cat per gene into the /dev/null sink), pchisq + data.frame 13%.
Bottlenecks: (1) large N: eigen set-up O(N^3) and rotation U'G O(N^2 S); (2) bootstrap ~45% at small N;
(3) Newton refine with GRM kernels (no eigenvalue ties); (4) per-gene R overhead (cat!) for many small genes;
(5) multi-core efficiency with few large genes; (6) RAM at large N (N x N copies in set-up).

## Round 12.1 — algorithm specialist (bench/proto/r12/algorithm/REPORT.md)
zz (analytic sparse Z Z' basis, default kernel): end-to-end 1.27x (S) .. 4.3x (XL) .. 21x (XXL), RSS -14..-57% at
N >= 1000; data l <= 1.2e-10; bootstrap draws change (needs distribution check, L KS p 0.010 -> second seed).
p12_schur (Schur reuse across models / bootstrap, identical()): 1.05-1.25x. p12_parab: methodological, 1.1x, held.

## Round 12.1b — reviewer audit and validation suite (bench/proto/r12/review/REPORT.md)
validate12.R (equiv / stat / repro / edge / stress; TSV + VERDICT): sanity base12 vs base12 PASS (1105 checks),
final vs base12 PASS; negative control MAX_PEAKS = 1 caught only by M_learn, edge and stress parts -> mandatory.
Audit fixes applied: P_many is really ~69 SNPs/gene (206,669 tests, overlapping windows; kept as a many-dense-genes
workload); new P_many10 (10 SNPs/gene), XL_grm_many (N = 2000 GRM, 400 genes), R_real (GRM + 10 covariates +
missing + uneven windows; ~20x M_cov cost per test) linked into bench/data; harness default -> base12; MKL /
FlexiBLAS / GOTO threads pinned in bench.sh; round12_bench.sh alternates implementation order per repetition.
Open risks for candidates: R1 unblocked per-gene O(N S) memory (883 MB for 4000 SNPs at N = 2000), R4 Newton on
GRM + covariates, R2 latent use-after-free in fit_models.
Baseline on the reviewer datasets (1 rep, machine contended by agents; indicative): R_real 198 s at 1 core (59,046
tests, 3.4 ms/test), 49.5 s at 4 cores, 352 MB; XL_grm_many 47.5 s / 516 MB; P_many10 7.5 s for 29,986 tests
(3000 genes: ~2.5 ms fixed cost per gene).

## Round 12.1c — parallel specialist (bench/proto/r12/parallel/REPORT.md)
P1 opt-in suez.threads (std::thread over SNP blocks; normals drawn on the main thread; identical() at 1-8 threads):
kernel 3.3-3.4x at 8 threads for XL_grm / XXL, ~1.0 at small N; end-to-end equal-ish to 8 forks at 28-31% of the
tree RAM. P2 macOS Accelerate dsyevd for explicit kernels: eigen 7x at N = 2000; XL_grm 1.53x at 1 core (draws
change). Rejected: Eigen+OpenMP (no libomp; <= 8% parallel), dynamic / largest-first scheduling (neutral).
zz validation (validate12.R algorithm_zz --boot dist, all 25 datasets): 966 PASS, 1 WARN (L bootstrap LR KS p 0.0096,
seed 42), 0 FAIL; data -log10 p <= 1.6e-8, lead SNPs / significant sets / eGenes identical, edge + stress PASS.
Three new seeds on L (bench/proto/r12/review/boot_seeds.R): KS p 0.061 / 0.35 / 0.67, Fisher pooled 0.21, quantiles
match to q99.9 -> the warning was chance. zz ACCEPTED (documented bootstrap-draw change for sample_kernel = NULL).

## Round 12.1d — numerical specialist (bench/proto/r12/numerical/REPORT.md)
p_prof2 (pkg_comb, lib numerical_comb, Eigen only, ~+150 lines): Profile as 3 GEMVs + faster LDLT with Eigen's
pivot order; refine reuses the best Newton call; grid covariate-block solve as one GEMM per column per SNP block;
vectorised moments; cheaper log-det on GRM; residual derivatives when rows ~ n (reviewer R4: the 3-column moment
product over all pairs on R ~ n rows was the cost). Kernel 1.13-1.17 (S), 1.2 (M), 1.4 (M_eqtl), 1.5 (M_cov),
1.8-2.6 (M_cov30), 1.2-1.5 (GRM), 1.6-1.9 (R_real), 1.04-1.18 (L/XL/XXL). Data l <= 1.6e-10, beta 6.3e-8, se 2.3e-8,
boot <= 8.4e-7 (same draws); LA 1.7e-13. validate12 quick: PASS with 3 warnings inside base12's own floor (one exactly
collinear SNP finite vs NA; near-perfect fit 8e-7; stress 1e-7, candidate closer to exact). Newton now starts from
Profile values (robustness fix: an extreme case returned the grid point 3.8e-9 below the max in both builds).
GEMM GFLOP/s one M5 core: Eigen 17-19, OpenBLAS 20-66, Accelerate 108-141, R reference BLAS 1.2 -> BLAS dispatch
only for optimised BLAS (Linux opt-in; macOS Accelerate a separate decision). Rejected: in-place moment columns,
residual derivatives on default kernels, unpivoted Cholesky, low-rank grid weights.
Next: c2 = c1 (zz + schur) + p_prof2, merged by the numerical specialist.
c1 = zz + p12_schur (bench/proto/r12/integ/pkg_c1, lib integ_c1): validate12 vs algorithm_zz --boot exact, all
datasets: PASS (960 PASS, 0 WARN/FAIL); every output column identical on every dataset (bootstrap included).

## Round 12.1e — memory specialist (bench/proto/r12/memory/REPORT.md)
All prototypes identical() to base12. Lazy per-chunk rotation in fit_gene: 829 -> 336 MB on a 4000-SNP gene at
N = 2000. Lean set-up .Call (pkg_exact, identical under any BLAS): XL_grm 0.76x, XXL 0.78x peak footprint.
R side: dplyr/readr not loaded in step 2 (library 124 -> 88 MB; tree RSS -12..-22% at 2-8 cores), cheaper per-gene
table. Combined memory_gene: XL_grm_many 0.54x, XL_grm 0.62x, M 0.95x footprint, wall 0.90-1.00x.
CORRECTION to round 12.0: cat() is not 38% of P_many (64 ms / 1000 calls measured directly); Rprof attributed the
next .Call's time to it.

## Round 12.2 — integration c2 = c1 + p_prof2 (numerical specialist; bench/proto/r12/integ/pkg_c2, lib integ_c2)
p_prof2 engine + c1's zz sparse path + bootstrap y-free moment reuse; c1's SchurCache folded into p_prof2's per-column
GEMM grid (data grid stores its y-free part per grid point, bootstrap grid reuses it); c1's final-value skip superseded
by p_prof2's best-Newton-call reuse. c2 identical() to numerical_comb on M_grm, L_grm, R_real.
validate12 vs integ_c1 --ds quick --boot exact: PASS with 3 warnings (same as p_prof2 vs base12, inside base12's floor);
data l <= 1.7e-10, beta 1.2e-7, boot 9.4e-7 (same draws).
Kernel (interleaved, 3 reps; base12 / c1 / c2 s): S 0.097/0.074/0.070, M 0.233/0.149/0.136, M_grm 0.529/0.502/0.400,
M_cov30 1.89/1.22/0.80, R_real (8 genes) 1.43/1.35/0.79, XL 0.97/0.25/0.23 -> base12/c2 1.32-4.19x, c1/c2 1.06-1.71x.
e2e 1 core: S 0.69/0.55/0.55 s, M_grm 3.38/3.10/2.42 s; RSS unchanged.
c2 full validation vs base12 (--ds all --boot dist): PASS with 4 known warnings (L KS seed 42 = zz, cleared by 3 seeds;
collinear edge case: an exactly collinear SNP finite vs NA, as base12 under a 2e-15 perturbation of y; near-perfect fit
7.95e-7 inside base12's own floor 4.4-7.6e-7; stress |c2 - base12| 1.06e-7 with c2 closer to the exact ML, shortfall
1.42e-13 vs 1.71e-13). Reproducibility identical() for rerun / 2 / 4 cores on S, M_dense, M_grm, LA.

## Round 12.3 — integration c3 = c2 + memory changes (memory specialist; pkg_c3, lib integ_c3)
Lazy per-chunk rotation (LazyRotated, same chunks as rotate_all -> bit-identical; fit_gene and fit_models; dense and
sparse U); explicit-kernel set-up: K H in R, then .Call(suez_contrast_eigen) running the eigen solver in place in U
(LAPACK path dsyevr as eigen()), remaining products through R's own dgemv -> identical under any BLAS; Makevars
PKG_LIBS $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS); R side: positional lookups, list-built per-gene table + bind_tables(),
dplyr/readr by namespace (not loaded in step 2). validate12 vs integ_c2 --ds all --boot exact: PASS, 0 WARN/FAIL;
identical() to c2 on all 41 cached outputs (28 datasets, null/none/permute, 2/4 cores).
Peak footprint base12 / c2 / c3 (MB, 1 core, median of 3): M 190/208/197, R_real 380/380/321, XL_grm 641/646/444,
XL_grm_many 642/644/460, XXL 1467/414/313. mem_gene.R (N = 2000): set-up 424/424/281, one 4000-SNP gene 809/828/371.
BLOCKER found: tests/model_checks.R kernel x1e6 / x1e-6 invariance FAILS in numerical_comb, c2, c3 (max |dl| 1.2e-5;
base12, zz, c1 pass with 8.6e-8) -> introduced by p_prof2; the validate12 suite did not cover it. Sent back to the
numerical specialist for diagnosis and fix (c3b). Lesson: run tests/*.R on every integrated candidate, not only validate12.
c3b (pkg_c3b, lib integ_c3b): regression cause = p_prof2's Newton start at -Inf instead of the grid value: a peak at
the grid boundary (h = 0 or 1) cannot be reached by Newton (x capped at +-12), so at kernel x1e6 a gene with optimum
h = 0 was fitted at delta = 1e-6 (l 3.7e-6..9.7e-6 below the exact ML on 8 of 48 rows; base12 exact to 4e-14).
Fix (1 line): boundary grid points start from P(m, hb). tests/*.R all PASS (model_checks 17/17, kernel-scale |dl|
8.6e-8 = base12); validate12 vs integ_c3 quick: PASS, 2 known warnings; cost <= 0.5%.
Build portability: $(FLIBS) removed from PKG_LIBS (Makevars, Makevars.win): it broke source installs on a Mac without
gfortran; LAPACK_LIBS / BLAS_LIBS are R's shared libraries carrying their own Fortran runtime. Plain R CMD INSTALL OK.
Suite extended (reviewer): kscale12.R (kernel scale + exact boundary optima) and a `tests` part; integ_c3 and
numerical_comb FAIL it, base12 and integ_c3b PASS. integ_c3b vs base12 full run: PASS with 3 known warnings, 0 FAIL
(1118 checks).

## Round 12.4 — c4 = c3b + opt-in threads (parallel specialist; pkg_c4, lib integ_c4)
Default identical() to c3b; identical() at 1/2/4/8 threads (S, M_grm, LA, R_real, XL_grm, XXL); validate12 vs c3b
PASS 0 WARN; tests PASS incl. a new thread-identity test. 1 process x 8 threads = 0.51-0.73x of 8 forks' speed at
20-24% of their tree RAM (XL_grm 516 vs 2183 MB); no gain on M. Integrated as opt-in (memory saver); forks default.
Large-N explicit-kernel walls barely moved (XL_grm 11.6 s vs 12.1 base12): set-up eigen + O(N^2 S) rotation on the
reference BLAS remain -> round 13 target.

## Round 12.5 — final serial benchmark and integration (coordinator)
bench/round12_bench.sh 3 base12 integ_c4 (serial, quiet machine, order alternated per rep;
bench/results/runs_r12_final.tsv, summary bench/results/r12_final_summary.tsv). Step 2, 1 core, wall base12 -> c4
(speed-up [min-max of reps]; peak RSS): S 0.66 -> 0.48 s (1.38x; -7%), M 2.92 -> 1.74 (1.68x; -6%), L 7.81 -> 3.15
(2.48x; -33%), XL 11.07 -> 2.23 (4.97x; -59%), XXL 30.56 -> 1.20 (25.4x; -78%), M_grm 1.40x, L_grm 1.36x (-19%),
XL_grm 12.18 -> 10.11 (1.20x; -30%), XL_grm_many 1.12x (-23%), M_cov 2.05x, M_cov30 2.29x, M_fac3 2.19x (-18%),
M_dense 1.64x (-27%), M_miss 1.26x, M_learn 1.19x, M_lc 1.22x, M_num 1.39x, M_eqtl 1.37x, M_unbal 1.68x,
S_rare 1.37x, HS 1.38x, LA 1.39x, P_many 1.42x, P_many10 1.56x, R_real 55.97 -> 35.07 (1.60x; -15%).
Multi-core: M 1.56-1.65x (tree RSS -17..-25%), M_dense 1.62-1.69x (-22..-34%), R_real 1.50-1.57x (-8..-16%),
XL_grm 1.06-1.16x (-7..-23%). Step 1 (code unchanged): 0.95x / 0.98x (noise). Minimum per-rep speed-up > 1 on every
step-2 row. Integrated pkg_c4 into src/, R/, tests/ (NAMESPACE, Makevars); repo install: tests/*.R 0 failures.
Reviewer datasets (P_many10, R_real, XL_grm_many) are regenerated by bench/proto/r12/review/datasets12.R.
Remaining large-N cost with explicit kernels (XL_grm): serial eigen-decomposition + O(N^2 S) rotation on the
reference BLAS (BLAS dispatch / Accelerate measured 2-7x on those parts but platform-specific: not integrated).
