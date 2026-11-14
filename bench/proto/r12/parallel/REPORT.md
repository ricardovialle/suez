# Round 12 — parallel specialist (transcribed by the coordinator from the hand-back)

Noise: load 10-21 during all timings (absolute walls 2-3x quiet: XL_grm 27.8 vs 12.1 s, XXL 85 vs 29.6 s);
interleaved, 3 reps, medians; ratios +-20%. Libraries: bench/lib/parallel_thr, parallel_acc.

## Where time goes at large N (1 core, base12; phases.R)
XL_grm (N 2000, 40 genes): set-up eigen 2.3 s, fits 9.6 s. XXL (N 4000, 16 genes): eigen 24.8 s, fits 12.4 s.

## Parallel set-up
Eigen + OpenMP: not buildable (no libomp); Eigen parallelises only Q formation (~8% of the solver; tridiagonalisation
28%, Givens QR 63%, eig_split.cpp) -> rejected. Linux already uses threaded LAPACK (round 11).
macOS Accelerate dsyevd on the Householder-projected kernel (eig_accel.*): N 1000 GRM Eigen 0.46 / dsyevd 0.07 s;
N 2000 GRM 5.12 / 0.73; N 2000 default 4.27 / 0.72; N 4000 ~24 / 14.0 (GRM), 8.5 (default). Eigenvalues 1.3e-12,
leading projector 7e-15. VECLIB threads made no difference: the gain is the LAPACK, not threads.

## Prototypes
P1 pkg_thr (thr.diff: lmm.cpp +113/-52, R +12/-2): options(suez.threads = n) (default 1); new .Call
suez_fit_gene_threads (8 args), suez_fit_gene unchanged. std::thread par_for (atomic counter), per-worker scratch /
Profile, exceptions re-thrown on main thread, no R API in workers. Rotation in 240-row chunks of the result
(Eigen GEMM bitwise identical for row splits in multiples of 24 on all shapes tested here; 64/256 rows or column
splits are not). Tasks = unchanged SNP blocks: data fits in parallel; normals drawn on the main thread in SNP and
contrast order; bootstrap refits in parallel. identical() to base12 at 1/2/4/8 threads on S, M, M_eqtl, M_fac3,
M_grm, XL_grm, XL, XXL and end-to-end for every forks x threads combination.
Kernel speed-up 2/4/8 threads: XL_grm 1.77/2.67/3.28, XXL 1.97/2.74/3.44, XL 1.52/2.00/2.16, M_grm ~1.2, small N ~1.0.
P2 pkg_acc = P1 + Accelerate dsyevd ("auto" on Darwin with reference LAPACK; options(suez.eigen_solver =
"accelerate")): data l <= 1.2e-10, beta 8.5e-8 rel, NA/lead SNPs/eGenes identical; bootstrap draws per seed change.
P3 scheduling (dynamic, largest-first) inside pkg_thr: identical outputs.

## End-to-end (contended; speed-up vs base12 1 core)
| ds | config | CPUs | wall | speed-up | tree RSS MB | identical |
|---|---|---|---|---|---|---|
| XL_grm | base12 1 core | 1 | 27.8 | 1.00 | 598 | - |
| XL_grm | base12 8 forks | 8 | 10.8 | 2.57 | 2579 | yes |
| XL_grm | P1 1x8 threads | 8 | 12.8 | 2.17 | 713 | yes |
| XL_grm | P1 4x2 | 8 | 9.6 | 2.90 | 1645 | yes |
| XL_grm | dynamic / lpt 8 forks | 8 | 12.9 / 11.2 | 2.16 / 2.47 | ~1800 | yes |
| XL_grm | P2 1x1 | 1 | 18.1 | 1.53 | 727 | data only |
| XL_grm | P2 1x8 | 8 | 7.2 | 3.85 | 775 | data only |
| XL_grm | P2 8 forks lpt | 8 | 4.7 | 5.94 | 1921 | data only |
| XXL | base12 1 core | 1 | 84.9 | 1.00 | 1042 | - |
| XXL | base12 8 forks | 8 | 59.5 | 1.43 | 3983 | yes |
| XXL | P1 1x8 | 8 | 67.8 | 1.25 | 1221 | yes |
| XXL | P2 1x1 | 1 | 34.7 | 2.45 | 1127 | data only |
| XXL | P2 1x8 | 8 | 20.9 | 4.07 | 1313 | data only |
| M | base12 8 forks | 8 | 2.36 | 3.85 | 1255 | yes |
| M | P1 1x8 | 8 | 8.64 | 1.05 | 246 | yes |

## Findings / recommendation
Threads win per byte, not per second (XL_grm 1x8 = 84% of 8 forks' speed at 28% RAM; XXL equal speed at 31% RAM;
M no gain). The large-N cap is the serial set-up. Scheduling changes neutral (uniform synthetic genes).
Adopt P2 after bootstrap-distribution validation or opt-in; P1 as opt-in suez.threads with a cross-thread-count
identical() test; drop OpenMP and dynamic/LPT scheduling. Not done: x86/GCC bit-invariance, R CMD check, fit_models
threading, quiet-machine timings.

## c4 = integ_c3b + P1 (opt-in threads) — bench/proto/r12/integ/pkg_c4, lib integ_c4 (diff: parallel/c4_vs_c3b.diff)
options(suez.threads = n) (env SUEZ_THREADS fallback, default 1); 8-argument suez_fit_gene_threads only when
threads > 1 (suez_fit_gene unchanged). SNP blocks in waves of nt: rotation (LazyRotated) on the main thread, data fits
in parallel (one Profile per worker), bootstrap Yb drawn on the main thread in c3b's order, refits in parallel. Dense U
rotation threaded over 240-row chunks inside the 256-column chunks (bit-identical: Eigen GEMM rows identical for
splits in multiples of 24); sparse rotation and fit_models serial. New test (model_checks.R): identical() at 1/2/4/8
threads with bootstrap, GRM (dense U, N = 600) and default kernel (sparse U).
Checks: tests 0 failures; validate12 vs integ_c3b --ds all --boot exact: PASS (973 PASS, 0 WARN/FAIL);
thread identity on S, M_grm, LA, R_real, XL_grm, XXL: identical().
Timing (3 interleaved reps, load 3-6; wall / tree RSS): M 1 core 2.08 s / 214 MB, 8 forks 0.52 / 1167, 1x8 threads
2.12 / 214, 4x2 0.76 / 702; XL_grm 11.6 / 478, 4.51 / 2183, 6.18 / 516, 5.08 / 1623; XL_grm_many 13.3 / 496,
5.45 / 2392, 10.0 / 495, 6.31 / 1433; R_real 38.9 / 381, 9.61 / 2339, 18.9 / 457, 11.1 / 1404.
=> 1x8 threads: 0.51-0.73x of 8 forks' speed at 20-24% of their RAM; forks stay the default; suez.threads documented
as a memory saver. Not checked: bit-identity of the 24-row split on x86/GCC (the new test guards it), R CMD check.
