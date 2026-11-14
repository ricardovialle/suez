# Round 12 — memory specialist (transcribed by the coordinator from the hand-back)

Caveat: machine heavily loaded (load 13-16, 4.9 GB swap): walls 2-3x baseline and noisy; RSS deflated by memory
compression, so RAM = peak phys_footprint (/usr/bin/time "peak memory footprint"). Medians of 3 paired interleaved
reps, 1 core unless stated. Files: pkg_* (prototypes), pkg_*.diff, scripts, runs_{eig,b,c,d,e}.tsv, setup_mem.tsv, out/.

## Ranked prototypes (all identical() to base12, bootstrap columns included)
1. memory_gene = pkg_gene (all below combined) — recommended. Peak footprint XL_grm_many 0.54x (627 -> 341 MB), XL_grm
   0.62x, XL_grm 8 cores 0.60x (parent; tree RSS 0.91x), M 0.95x, P_many10 1.01x; wall 0.90-1.00x. validate12
   memory_gene --ds quick PASS (0 FAIL); tests 35 PASS on memory_all3 / memory_exact2 (not re-run on memory_gene).
2. Gene-loop lazy rotation (reviewer R1; C++ only, ~40 lines in fit_gene): G not copied (Map; integer G accepted);
   Gs / Is standardised and rotated per 256-column chunk with exactly rotate_all's chunks (bit-identical), <= 2 chunks
   held. mem_gene.R (N = 2000, S = 4000): peak 829 -> 336 MB; per-gene working set +405 -> +85 MB; fit 2.17 vs 2.02 s.
   Not applied to fit_models yet.
3. Lean set-up (pkg_eig / pkg_all): one .Call(suez_contrast_eigen) forms H K H from K on the fly (reference-dgemv sum
   order with FMA), Eigen's SelfAdjointEigenSolver::compute replayed in place in the returned U; LAPACK path calls
   dsyevr as La_rs with ldz = N. Extra memory: Eigen path 6.9 -> 1.9 N^2 doubles, LAPACK path 4.9-6.9 -> 2.9 N^2.
   End-to-end XXL 0.63x (1471 -> 921 MB), XL 0.66x, XL_grm 0.70x. Identical only under the reference BLAS (OpenBLAS /
   MKL: eigenbasis differs at ulp level -> bootstrap draws per seed change). pkg_exact: H K H and u' products through
   R's BLAS calls -> identical under any BLAS for +1 N^2 (2.9 / 3.9 N^2; XL_grm 0.76x, XXL 0.78x). Recommended if
   Linux bit-identity matters.
4. R side (pkg_rside + bind_tables): gene coordinates by position, double round, exact constant-column test, per-gene
   table built as a list, own row binding; dplyr/readr called by namespace and not loaded in step 2: library(suez)
   124 -> 88 MB; multi-core tree RSS -12..-22% (M 2/4/8 cores), wall 0.94-0.97x; per-gene R overhead -0.9 ms
   (P_many 0.85-0.98x).

## Findings
- cat() is NOT 38% of P_many: 64 ms per 1000 calls when timed directly; Rprof attributed the following .Call to it.
  data.frame() 0.5 ms/gene; pchisq cheap.
- Forked workers: U not copied (header page only, read-only Map); workers copy small-object heap pages dirtied by GC
  -> fewer loaded namespaces = less copied.
- Result table ~150 B/row of per-gene lists + ~130 B/row final; column-wise freeing does not lower the peak; only
  streaming would.
- Rejected: row-blocked UZ (not bit-identical, 1.3e-15), gc() per gene, avoiding Z.
Not checked: GCC/Linux (no gfortran here; PKG_LIBS uses $(LAPACK_LIBS) $(BLAS_LIBS) without $(FLIBS)), bootstrap
distribution under OpenBLAS, the exact path on a real optimised BLAS.
