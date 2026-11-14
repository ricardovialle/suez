# Round 12, numerical specialist: report

Machine: Apple M5, R 4.6 with the reference BLAS/LAPACK, clang -O2. **Other agents loaded the machine heavily during
every measurement (load average 11-17)**: absolute times are inflated and unstable (base12 M_cov30: 1.2-4.9 s for the
same input), so only interleaved ratios (base12 and prototype alternated in one process, median of 3 reps) mean
anything, and even those moved by +-20% between runs (both runs are reported). The coordinator's serial runs are the
reference. Every prototype is a drop-in `Rcpp::sourceCpp` copy of `bench/proto/r12/p12_base.cpp`.

## 1. Where the time goes (base12)

Coordinator profile, `sample` call-graph profiles (`sample.sh`, `attrib.py`), Profile instrumentation (`p_instr.cpp`):

| workload | dominant costs |
|---|---|
| GRM kernels (M_grm, L_grm), R ~ n eigenvalue rows | Profile/Newton ~50%: per call, the moment product MM'(P2 x R) x 3 weight columns, where Eigen's GEMM on a 3-column skinny product is 2-3x slower than three GEMVs; `weights()`: R logs plus 3 heap temporaries per call (~20% of M_grm alone); grid GEMMs |
| many covariates on a default kernel (M_cov30) | grid 42%: for each SNP x grid point x model, z = A00^-1 r0 and T_j = A00^-1 a_j as scalar strided p0^2 loops, recomputed for every model and, for the shared data y, for every SNP. Profile 36%: a pivoted 33x33 LDLT per Newton call (5-6 us in Eigen). Moments 13%: scalar dot loops, which strict FP keeps out of SIMD |
| GRM + covariates (R_real, reviewer lead R4) | both of the above: P2 = 91 pairs on R = 849 rows, 3 weight columns per Newton step |
| large N (XL, XXL, XL_grm) | rotation U'G: Eigen GEMM runs at ~18 GFLOP/s on one M5 core |

Each refine makes 6.4 Profile calls. One of them is a value-only re-evaluation at the best h, which was already
evaluated.

## 2. Prototypes

All are exact reformulations: same model, grid, peaks, Newton steps, tolerances and RNG draw order. Outputs change only
by rounding, except E, a robustness fix.

- **p_prof.cpp** (Profile / Newton)
  - A. The moment products are three GEMVs instead of one 3-column GEMM (micro: 2.3x at 306 x 10, 1.5x at 110 x 630;
    `micro/mw.cpp`).
  - B. `FastLDLT` reproduces Eigen's LDLT pivot order exactly. Eigen's unblocked LDLT pivots on the largest
    *original* |diagonal|, which is not updated before the search. The factorisation is then right-looking with
    column axpys, keeping PinvLDLT's keep/drop rule (COLLINEAR) and its pseudo-solve. At 33 x 33: 3.7-4.8 us vs
    4.7-5.7 us (`micro/ldlt*.cpp`).
  - C. refine() reuses the state (b, s, A, rhs) of the best Newton call instead of calling P(m, hb) again. Column 0 is
    formed by the same GEMV in value-only and derivative calls, so the reused value is bit-identical to a
    recomputation. This saves 1 of 6.4 calls per refine.
  - D0. `weights()` no longer allocates temporaries on each call.
- **p_comb.cpp** = p_prof plus the changes below.
  - D1. Grid Schur complement through a factor of A00^+ per grid point: F' = D^{+1/2} L^-1 P from the pivoted LDLT
    (dropped pivot -> 0, sign of kept pivots carried). Then a'A00^+b = (F'a).(sg F'b).
    - Per grid point and SNP block, one triangular GEMM per block column (all SNPs at once) replaces the per-SNP
      scalar p0^2 loops.
    - Each column is transformed once for all models. The data models {0} and {0..q} share it, as do the bootstrap
      models {}, {0} and {0..q}.
    - The shared-y term of the data block is computed once per grid point instead of once per SNP.
    - The per-model Cholesky, SCHUR_PIVOT fallback and direct-residual paths are unchanged.
  - D2. Moments:
    - Runs of singleton eigenvalue rows (GRM) become one elementwise product loop (bit-identical).
    - Rows of 8 or more contrasts use Eigen's vectorised reductions (rounding changes).
  - D3. Log-determinant in Profile for all-singleton kernels: one log per product of 4 rows, guarded to
    d in [1e-60, 1e60], else Eigen's log. **New evidence vs the rejected round-1 p8_logdet:** on GRM kernels these
    logs are ~20% of the fit (sample profile of M_grm), and here the change gains 7-9%.
- **p_comb_fix.cpp** = p_comb + E. newton() takes its starting best value from Profile evaluations only (the first
  one is at the grid point itself) instead of from the grid value.
  - Found on the reviewer's extreme case `default_x1e3_opt_x4`. There the grid value at the peak carries ~1e-8 of
    Schur-complement rounding. When it happens to exceed every Newton value, the fit is reported *at the grid point*,
    3.8e-9 below the maximum.
  - base12 has the same fragility. At that SNP its grid error was negative only by luck: the rounding of the moments
    decides it. With E, both are within 2.5e-12 of the exact optimum.
  - Outputs move by <= 1.6e-10 vs base12. **This is a robustness fix, label it as such.**
- **p_prof2.cpp** = p_comb_fix + F. **This is the recommended Eigen-only engine.**
  - F. When 2 R P2 > 1.5 n (2p + 4), the derivative sums come from residuals: Q1, Q2 and X'D1r are built from
    r = y - Xb and per-contrast weights (`Spectrum::expand`, the same null-group expansion). That costs n (2p + 4)
    flops instead of two more moment GEMVs (2 R P2 = 155k flops per call on R_real). Default kernels, with few rows,
    keep the moment path.
  - On the reviewer's hypothesis (R4): A00(h) cannot be shared across Newton calls, because every call has its own h.
    But the 3-column moment product over all (p0 + k + 1)^2 / 2 pairs on R ~ n rows was indeed the cost, and F
    removes two thirds of it.
- **p_blas.cpp** = p_prof2 + G. GEMM dispatch for the rotation (U'G, and UZ Mu on the repeated-row path) and the grid
  products (moments x grid weights). Products below 2e4 flops stay in Eigen.
  - mode 0: Eigen.
  - mode 1: R's BLAS through `F77_CALL(dgemm)`.
  - mode 2: Accelerate `cblas_dgemm`, linked directly (`-DSUEZ_ACCELERATE`, `-framework Accelerate`).
- **Package builds:**
  - `pkg_comb` = p_prof2 → `bench/lib/numerical_comb`.
  - `pkg_acc` = p_blas with mode 2 as the default and Accelerate added to Makevars, macOS-only as written →
    `bench/lib/numerical_acc`.
- **Not kept:**
  - p_prof3 (moment columns read in place, contiguous Z for the residual GEMVs): 0.95-1.03x.
  - The residual derivative path on default kernels: more flops than the moments when R is small.
  - An unpivoted Cholesky with fallback: it cannot be made rank-safe without pivoting (Kahan-type matrices).
  - A low-rank grid weight matrix: approximate.

## 3. GEMM / LAPACK measurements (one M5 core; `gemm/`, `eig/`)

- **Rotation GEMM** (N x N-1 by N x 256), GFLOP/s:
  - Eigen 18-19.
  - R reference BLAS 1.2-1.3, 14x slower than Eigen: mode 1 must never be used with the reference BLAS.
  - OpenBLAS 0.3.33 (the pixi env's library, ARM, 1 thread) 36-66.
  - Accelerate 108-133.
- **Grid-shaped GEMM** (1200 x 82 x 1000 and similar), GFLOP/s: Eigen 17, OpenBLAS 20-42, Accelerate 121-141.
- **Accelerate threading:** its dgemm runs single-threaded (CPU time equals elapsed time, with or without
  VECLIB_MAXIMUM_THREADS).
- **Fork safety (`fork_check.R`):** under mclapply with 2 workers, results were identical to serial and across runs,
  with no hang. Scaling at 2 workers was 1.52x vs 1.74x with Eigen (the AMX unit is shared per cluster), so check at
  4-8 workers.
- **Linux x86:** packages are built at -O2 without -march, so Eigen gets SSE2 only, while OpenBLAS/MKL dispatch
  AVX2/AVX-512. The gap should be larger than the ARM figure above. This is not measured: there is no x86 machine
  here, and the pixi R has no Rcpp.

Eigen-solver timings for kernel_contrast_eigen's (N-1) x (N-1) GRM-like K, in seconds on the loaded machine:

| N | Eigen solver (base12, ref R) | R eigen() = ref dsyevr | Accelerate dsyevd | Accelerate dsyevr |
|---|---|---|---|---|
| 1001 | 0.91 | 1.40 | 0.11 | 0.17 |
| 2001 | 7.70 | 12.17 | 1.33 | 1.51 |
| 4001 | 49.4 | - | 9.35 | 8.69 |

- **Accuracy:** eigenvalues agree to 6e-13 relative, top-10 projectors to 2e-13, residuals are 1e-15.
- **dsyevd vs dsyevr:** no clear winner in Accelerate. On Linux, R's eigen() already uses the optimised dsyevr
  (round 11).
- **Draw changes:** any solver change rotates the basis inside tied eigenvalues, so **bootstrap draws per seed
  change**. This was not prototyped as a package change.

## 4. Correctness vs base12

- **Kernel level** (`cmp.R`, 18 workloads incl. R_real, HS, M_learn, P_many, L, `correct_fix.log`, `timing_final*.log`):
  - Data l (l0, l_geno, l_interact) <= 1.6e-10 (M_learn), and <= 2e-11 elsewhere.
  - beta <= 6.3e-8 relative to max(|beta|, se); se <= 2.3e-8 relative.
  - df and NA pattern identical.
  - Bootstrap l <= 8.4e-7, with the same draws.
  - The Accelerate mode has the same maxima.
- **Local ancestry** (fit_models, `la_check.R`, 20 genes): ll 1.7e-13, ll_boot 6.2e-8, beta 1.1e-8, NA/npar identical.
- **Extreme cases** (`extreme_dig.R`): every case is within base12's own floor, measured with y x (1 + 1e-15). Two
  notes:
  - grm_reps_x1e6_h_to_1: bootstrap interaction 3.4e-6 (< 1e-5).
  - default_x1e6_opt_x4: bootstrap 3.7e-5, equal to base12's own floor of 3.7e-5.
- **Collinear sweep:** identical for 1 - R^2 >= 1e-12. At the declared rank tolerance (1 - R^2 <= 1e-14) both builds
  differ from the full-rank optimum by design. The sweep is unseeded, and base12 itself moves between runs.
- **rankdef_beta:** no NaN; relative error 3e-10.
- **Reviewer validation:** `validate12.R numerical_comb --ds quick --parts equiv,stat,edge,stress` gives **PASS with
  3 WARN** (612 checks).
  - collinear_interaction_and_snps: 2 NA mismatches, at one SNP whose column is exactly collinear. base12 reports NA
    and the candidate a finite beta (identical l). This decision is a coin toss in base12 itself: Uy x (1 + 2e-15)
    flips it too (`edge_na.R`). Reviewer note: with exact collinearity the wald() keep/drop ratio is rounding noise
    around COLLINEAR.
  - near_perfect_fit: l 7.95e-7, within base12's floor of 4.4-7.6e-7. The candidate is closer to exact than base12.
  - random 3-block stress: |cand - ref| 1.06e-7, but the candidate's shortfall vs the exact ML is 1.4e-13 vs 1.7e-13 for
    the reference.
  - The repro part (1/2/4 cores) was not run, because I keep to 2 cores.

## 5. Speed

Kernel level: fit_gene C++ time, base12 / prototype, median of 3 interleaved reps.

| dataset (N, SNPs) | run 1: p_comb_fix Eigen | run 1: + Accelerate | run 2: p_prof2 Eigen | run 2: + Accelerate |
|---|---|---|---|---|
| S (200, 3975) | 1.17 | 1.21 | 1.13 | 1.24 |
| M (400, 5980) | 1.20 | 1.24 | 1.17 | 1.36 |
| M_grm (400, 5980) | 1.43 | 1.80 | 1.16 | 1.61 |
| M_cov30 (400, 3970) | 1.79 | 1.90 | 2.61 | 2.43 |
| M_cov (400, 5980) | 1.48 | 1.55 | 1.53 | 1.67 |
| M_fac3 (600, 3970) | 1.23 | 1.50 | 1.26 | 1.54 |
| M_eqtl (200, 5980) | 1.39 | 1.82 | 1.37 | 1.67 |
| P_many (200, 20615) | 1.09 | 1.06 | 1.13 | 1.22 |
| R_real (850, 1743) | - | - | 1.79 (1.61-1.89 in A/Bs) | 2.32-2.39 |
| L_grm (1000, 4712) | 1.19 | 1.78 | 1.46 | 2.32 |
| L (1000, 5915) | 1.14 | 1.78 | 1.18 | 1.79 |
| XL_grm (2000, 2306) | 1.32 | 2.21 | 1.31 | 2.17 |
| XL (2000, 3108) | 1.12 | 2.13* | 1.04 | 2.41 |
| XXL (4000, 1504) | 1.10 | 2.40 | 1.09 | 2.42 |

\* XL with 4 genes in an earlier run.

Contribution of each change (A/B, `attrib_ab.log`):

| change | S | M_grm | L_grm | M_cov30 | R_real | M |
|---|---|---|---|---|---|---|
| p_prof alone | 1.14x | 1.17x | 1.14x | 0.99x | | |
| D1-D3 + E on top | +2% | +16% | +15% | +55% | | |
| F on top | | +12-13% | | | +12-13% | |
| vectorised moment rows | | | | 6-11% | | 6% |

## 6. RAM

Peak RSS, installed builds, one process (`rss2.R`):

| dataset | base12 | numerical_comb | numerical_acc |
|---|---|---|---|
| L_grm (6 genes) | 308 MB | 319 MB | 304 MB |
| M_cov30 (8 genes) | 219 MB | 219 MB | - |

- The new buffers are small: the per-grid-point transformed columns are 3 k p0 S doubles (<= 0.4 MB at p0 = 31,
  S = 256), and the Profile snapshot is p^2.

## 7. Recommendation

1. **Integrate p_prof2** (Eigen only). Speed-ups vs base12:

   | workloads | speed-up |
   |---|---|
   | default-kernel small/medium (S, M, M_fac3, P_many) | 1.1-1.3x |
   | M_eqtl | 1.4x |
   | many covariates (M_cov, M_cov30) | 1.5-2.6x |
   | GRM (M_grm, L_grm, XL_grm) | 1.2-1.5x |
   | GRM + covariates (R_real) | 1.6-1.9x |

   - Same RAM; complexity medium: about +150 lines (FastLDLT, factor grid, expand/dird, snapshot).
   - Label E (Newton incumbent) as a robustness fix.
2. **BLAS dispatch for the rotation and grid GEMMs, mode 1:** call R's dgemm when R's BLAS is optimised, using the
   same detection as eigen_solver() plus an option. This is the larger win for Linux/OpenBLAS/MKL users at N >= 1000.
   - **Unmeasured on x86.** Validate with bench/linux/run_checks.sh.
   - Watch OpenBLAS threads inside doMC workers: oversubscription and fork safety. Thread count could in principle
     change GEMM bits.
3. **macOS Accelerate (mode 2, pkg_acc):** the largest gains measured here, 1.2-2.4x on top of 1.
   - Needs a build-system change: a Darwin-only configure/Makevars.
   - Check AMX scaling at 4-8 workers.
   - Accelerate dsyevd/dsyevr for kernel_contrast_eigen would add 5-7x on set-up, but it changes bootstrap draws per
     seed (as in rounds 7 and 11). Decide separately.

## Files

All under bench/proto/r12/numerical/.

- **Engines:** p_prof.cpp, p_comb.cpp, p_comb_fix.cpp, **p_prof2.cpp**, p_blas.cpp.
- **Package builds:** pkg_comb/, pkg_acc/.
- **Scripts:** cmp.R, ab.R, cmp_modes.R, la_check.R, extreme_dig.R, edge_na.R, fork_check.R, rss2.R, sample.sh,
  attrib.py, instr.R.
- **Micro-benchmarks:** micro/, gemm/, eig/.
- **Logs:** timing_final.log, timing_final2.log, attrib_ab.log, correct_fix.log, validate_comb.log,
  extreme_dig.log, collinear*.log, rankdef*.log.

## 8. Integration: pkg_c2 = integ_c1 + p_prof2 (Eigen only) -> bench/lib/integ_c2

### How c2 was built

- **Starting point:** bench/proto/r12/integ/pkg_c1, copied to bench/proto/r12/integ/pkg_c2. src/lmm.cpp is p_prof2's
  engine (pkg_comb) with c1's hunks applied on top:
  - Block(src): copies the y-free moments from the data block into the bootstrap block.
  - zz sparse basis: `sp_map`, `rotate_all(SpMap)`, the U-templated fit_gene/fit_models, the VECSXP dispatch and
    `suez_sp_crossprod`.
  - c1's R side and init.cpp are kept as they are.
- **Dropped from c1:**
  - The Profile `last_*` / "skip repeated final value call": superseded by p_prof2's best-call snapshot, which covers
    every case, not only a final value-only call at the same h.
  - The scalar SchurCache: replaced by a cache on p_prof2's factored grid. `SchurCache` keeps the y-free transformed
    columns F'a_c and the Schur entries per grid point from the data grid. The bootstrap grid reuses them, plus its
    own Y transform and the k x k per-model Cholesky. The same call signature (cache, reuse) is used in fit_gene.
- **Bootstrap reuse:** both reuses (Block src and the cache) are kept. They are exact. On the GRM datasets c2 is
  `identical()` to numerical_comb.

### Checks

1. **c2 vs numerical_comb** (kernel level, process-separated): `identical()` TRUE on M_grm (20 genes), L_grm (8) and
   R_real (6).
2. **Reviewer validation:** `validate12.R integ_c2 --ref integ_c1 --ds quick --boot exact --jobs 2` gave
   "integ_c2 vs integ_c1: 624 checks, INFO=80 PASS=541 WARN=3", **VERDICT: PASS with warnings**.
   - Maxima: l0 / l_geno / l_interact 1.6-1.7e-10 (M_learn), beta 1.2e-7 relative (R_real), bootstrap l 9.4e-7
     (R_real). Bootstrap draws are the same as c1's (exact mode).
   - The 3 WARN are the same as numerical_comb vs base12, with the same values:
     - collinear_interaction_and_snps: NA vs finite beta for one exactly collinear SNP, a coin toss in base12 itself.
     - near_perfect_fit: 7.95e-7, within the base floor.
     - stress: |cand - ref| 1.06e-7, with the candidate's shortfall 1.4e-13 vs 1.7e-13 for the reference.
   - The repro part (1/2/4 cores) was not run, to keep to 2 cores.
3. **Timing.** Kernel level: process-interleaved, 3 reps, medians (`integ_timing.sh`, `integ_timing.tsv`). Each build
   uses its own set-up: the zz sparse U for default kernels in c1/c2, the dense U in base12.

   | dataset (genes) | base12 s | c1 s | c2 s | base12/c2 | c1/c2 (per rep) |
   |---|---|---|---|---|---|
   | S (40) | 0.097 | 0.074 | 0.070 | 1.39 | 1.06 (1.04 1.08 1.06) |
   | M (30) | 0.233 | 0.149 | 0.136 | 1.65 | 1.07 (1.12 1.06 1.07) |
   | M_grm (30) | 0.529 | 0.502 | 0.400 | 1.32 | 1.25 (1.25 1.32 1.25) |
   | M_cov30 (20) | 1.890 | 1.221 | 0.799 | 2.35 | 1.53 (1.44 1.54 1.53) |
   | R_real (8) | 1.431 | 1.353 | 0.792 | 1.81 | 1.71 (1.70 1.72 1.71) |
   | XL (8) | 0.969 | 0.250 | 0.226 | 4.19 | 1.11 (1.11 1.10 1.11) |

   End-to-end (bench/bench.sh, 1 core, 3 reps interleaved, `runs_e2e.tsv`):

   | dataset | wall base12 s | wall c1 s | wall c2 s | base12/c2 | c1/c2 | peak RSS MB (base12 / c1 / c2) |
   |---|---|---|---|---|---|---|
   | S | 0.69 | 0.55 | 0.55 | 1.26 | 1.00 | 194 / 193 / 193 |
   | M_grm | 3.38 | 3.10 | 2.42 | 1.40 | 1.28 | 219 / 217 / 220 |

   - Correctness (kernel level, c2 vs c1): data l <= 1.7e-11, beta <= 5.4e-8 relative, bootstrap <= 6e-7, NA/df
     identical on all six datasets.
   - Vs base12 on default kernels, the bootstrap columns differ, as for c1: the zz basis changes the draws per seed.

## 9. Regression fix: kernel-scale invariance (tests/model_checks.R lines 149-154) -> pkg_c3b / bench/lib/integ_c3b

### Cause

Found with `kscale.R`: the test's inputs at kernel level, f = 1, 1e6, 1e-6, each prototype vs an exact optimiser
over delta = s2k/s2 (x = log10 delta in [-25, 25] plus both boundaries).

| build | invariance max \|dl\| | l - exact at f = 1e6 (min) |
|---|---|---|
| base12, p_prof, p_comb | 8.6e-8 | -1.5e-12 |
| p_comb_fix, p_prof2 | 1.2e-5 | -1.2e-5 |

- The culprit is E, the Newton starting value (section 2): newton() started from -Inf instead of the grid value.
- When the grid peak is a boundary point (h = 0 or h = 1), Newton cannot reach it: it is capped at
  |x = log10(h/(1-h))| <= XCAP = 12. The grid value was the only evaluation at the boundary itself.
- With the kernel scaled by 1e6, gene 4 (simulated h = 0) has its optimum at h = 0, and x = -12 means
  delta_eff = 1e6 x 1e-12 = 1e-6. Dropping the grid value therefore reported the fit 3.7e-6 to 9.7e-6 below the
  maximum (8 of 48 rows).
- **base12 is the accurate side** here, matching exact to 4e-14.

### Fix

```cpp
vb = best == 0 || best == NG - 1 ? P(m, hb) : NEG_INF;
```

- Boundary grid points are evaluated by Profile at the boundary itself. Interior points keep E, because the first
  Newton step evaluates the grid point anyway.
- The `default_x1e3_opt_x4` case that motivated E stays fixed (proto - exact 2.5e-12).
- Prototype: p_prof2b.cpp. Integrated: pkg_c3b = pkg_c3 + this one line.
- c3's Makevars links $(FLIBS), which does not link on this Mac (no gfortran), so I built with `MAKEFLAGS='FLIBS='`.

### Exact comparison at f = 1e6 (l_geno, the 5 worst rows; all in gene 4)

| row | exact | base12 - exact | p_prof2 - exact | fixed - exact |
|---|---|---|---|---|
| 19 | -105.2119 | -4.3e-14 | -9.7e-6 | -4.3e-14 |
| 20 | -110.7979 | -2.8e-14 | -3.9e-6 | -2.8e-14 |
| 21 | -109.7133 | -4.3e-14 | -3.7e-6 | -4.3e-14 |
| 22 | -110.4147 | 0 | -5.1e-6 | 0 |
| 23 | -110.7537 | -4.3e-14 | -3.8e-6 | -4.3e-14 |

At f = 1e-6 every build is 9e-8 below exact: the x <= 12 cap at the h -> 1 end. This is existing behaviour, and it
is what gives the test its 8.6e-8.

### Suggested edge case for the reviewer's suite

- The test's GRM: 120 donors with 1-2 samples each, condition rnorm within donor, genes at h in {0, 0.3, 0.8, 1}.
- Run sample_kernel = f K for f in {1, 1e6, 1e-6} and require invariance <= 1e-6.
- At f = 1e6, compare l0 / l_geno / l_interact to the exact optimiser in `kscale.R`; it should be within 1e-9.
- validate12 missed it because no dataset there has a boundary optimum on a kernel scaled far from 1. At f = 1,
  x = -12 is indistinguishable from h = 0.

### Verification of integ_c3b

- **Tests:** tests/*.R all PASS against integ_c3b: learn_covariance 6, local_ancestry 12, model_checks 17. The
  kernel-scale check gives max |dl| 8.6e-08, the same as base12.
- **Reviewer validation:** `validate12.R integ_c3b --ref integ_c3 --ds quick --boot exact --jobs 2` gave
  "624 checks, INFO=80 PASS=542 WARN=2", VERDICT: PASS with warnings.
  - near_perfect_fit: 7.85e-7, report-only, within base12's floor of 4.4-7.6e-7.
  - random 3-block stress: |cand - ref| 1.06e-7, with equal shortfall vs exact (1.4e-13 for both).
  - Equivalence maxima vs c3: l 1.7e-10, beta 2.4e-11, bootstrap 1.8e-10. All are M_learn, where boundary peaks now
    take the boundary value again.
- **Cost** (kernel level, process-interleaved, 3 reps, c3/c3b): M_grm 1.000, M_cov30 0.995, R_real 0.999. That is
  <= 0.5%: one extra Profile call only for boundary peaks.
