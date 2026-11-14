# Algorithm specialist — round 12 (report transcribed by the coordinator from the agent's hand-back)

Timings on a loaded machine (absolute times 2-4x the coordinator's); all ratios base12 / prototype, interleaved,
median of >= 3 reps. Nothing in src/, R/, tests/ edited.

## 1. zz — analytic sparse basis for the default kernel Z Z' (sample_kernel = NULL)  [recommended]
On the contrasts orthogonal to the intercept Z Z' has three invariant subspaces: within-individual contrasts
(eigenvalue exactly 0); between individuals with the same sample count s (eigenvalue exactly s); between size
classes (k - 1 eigenvalues of a k x k problem, k = number of distinct sample counts). A Haar-like orthonormal
basis of each tied subspace makes U sparse (~N log2 N nonzeros: 48k at N = 4000 vs 16M). Basis built in 0.02 s
at N = 4000 (base12 eigen 17 s in the coordinator profile). Rotation U'G is O(nnz S); exact zero null
eigenvalues shrink the null group from 4 moment rows to 1.
End-to-end 1 core: S 1.27x, M 1.34x, M_unbal 1.40x, L 2.16x, XL 4.32x, XXL 21.3x; max RSS L 378 -> 325, XL 670 ->
366, XXL 936 -> 404 MB (small N unchanged). Kernel level 1.22x (S) - 6.11x (XXL).
Data columns vs base12: l 5.7e-13 (S) .. 1.2e-10 (XXL); beta (relative to max(|beta|, se)) <= 1e-7, se <= 2e-8;
df, NA pattern, lead SNPs, eGenes identical. Bootstrap draws per seed change (basis inside tied eigenvalues; as
round 7): null LR KS p 0.58, 0.80, 0.18, 0.010 (L), 0.41, 0.56; quantiles match to q99.9 — L needs a second-seed
check. Package pkg_zz/ (~50 lines), lib bench/lib/algorithm_zz; tests model_checks + local_ancestry_checks 0 fail.
Applies only to sample_kernel = NULL (an explicit Z Z' is not detected).

## 2. p12_schur.cpp — exact, identical() on 6 datasets  [recommended]
y-free Schur factor of the largest model computed once per (h, SNP) and reused by the smaller model and the
bootstrap grid; bootstrap reuses the data block's y-free moments; refine skips a repeated final value call.
M_cov30 1.25x, M_cov 1.20x, M_fac3 1.18x, M_grm 1.06x, S 1.05x, M_eqtl 0.97x (noise). ~60 lines; fit_models not
patched. p12_boot.cpp alone (bootstrap moment reuse): 1.03-1.09x, identical.

## 3. p12_parab.cpp — METHODOLOGICAL (Newton starts at the 3-point parabola vertex)
Derivative calls per refine 4.4 -> 3.8-3.9; 1.10-1.15x on S, M_grm; data l <= 3.7e-13, boot l <= 7.7e-7 vs
base12. Changes the iteration path: separate decision only.

No other exact Newton savings (weights depend on h; each model visits different h). Low-rank grid weights not
prototyped (cancellation risk, <= ~10% of time).
Files: zz_basis.R, test_basis.R, p12_zz.cpp, p12_boot.cpp, p12_combo.cpp, p12_newton.cpp, *_count.cpp,
bench_zz.R, bench_cpp.R, e2e.sh, boot_check.R, zz_*.log, runs.tsv.
