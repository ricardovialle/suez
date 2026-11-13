# Algorithm specialist — report (saved by the coordinator from the agent's hand-back)

**Recommendation: p3_compress_fast.cpp** — same model (grid 80 + h=0,1; MAX_PEAKS 3; Newton |x|<=12;
ok_nested; RNG order), 2-7x faster at kernel level.

## p3: exact compression of the n contrasts
- Group contrasts by eigenvalue (relative 1e-12); null group (lambda <= 1e-12 max) handled with a 2nd-order
  Taylor expansion in delta = lambda - c_g, checked per h ((h max delta / d)^3 <= 1e-13), else falls back to the
  baseline per-contrast Profile. Singletons stay as they are.
- Per-group QR of [UX0, UG_s, I_s, y] keeps direct-residual arithmetic (no normal-equation cancellation).
- One compression per (SNP, y) serves all nested models (prefix columns); grid = one GEMM; one augmented
  Cholesky per (h, SNP) gives RSS of all nested models. Cancellation guard recomputes grid values with direct
  residuals when rss < 1e-10 y'Wy. Newton's final iterate reused. beta/se via the same LU inverse.
- Compressed rows m: default kernel 16 (n 199-999); M_unbal ~40; learned ~100; GRM ~n/2.

Correctness vs baseline (17 datasets): l0/l_geno/l_interact <= 1.6e-12; l_boot 7e-8 .. 8.6e-7 (L_grm);
identical NA/df; |dbeta|/max(|beta|,se) <= 5e-8. Baseline self-noise floor (Uy * (1 + 1e-15 N(0,1))): l_boot up
to 7.6e-7, dbeta/se 4e-8 — Newton's argmax over near-equal values defines hb only to ~1e-8, the bootstrap draws
inherit it. edge_cases.R: 21 PASS, 2 DIFF on report-only cases (duplicated covariate: baseline NaN, p3 finite
arbitrary; near-perfect fit: p3 matches an exact QR reference, baseline off by 6.6).

Speed (CPU ratio, interleaved vs sourceCpp p0, median 3):
S 5.6x, S_rare 5.6x, M_fac3 6.8x, L 6.5x, M_unbal 5.6x, M_learn 5.6x, M_lc 5.9x, M_cov 3.7x, M_grm 3.1x,
HS 3.0x, L_grm 3.1x, M_eqtl 1.8x. (Not interleaved: XL 7.6x, M_dense 7.2x, M_miss 6.5x, M_num 3.1x.)
Remaining time: rotation U'G now dominant (61% on L, 29% L_grm).

## Rejected
- Normal equations in Newton on grouped sums (unstable; QR compression makes it unnecessary).
- Snapping null eigenvalues to 0 (changes the model ~1e-4 in weights at x = 12).

## Householder kernel_contrast_eigen (eig_householder.R)
O(n^2) instead of O(n^3) products: saves 13.2 s on XL (21% of 62 s), 1.5 s on L. eigen() unchanged.
But a rounding-level change rotates the basis inside repeated eigenvalues: data fits identical (<= 2.6e-12),
bootstrap draws change completely (l_boot |diff| up to 11). Not a pure speed-up under the same-draws rule.

## Rotation caching across overlapping cis windows (overlap.R, not prototyped)
Mean overlap 2.0 genes per SNP (9.9 on M_dense). Est. saving after p3: ~30% L, ~15% L_grm, ~50% M_dense.
Requires gene chunks in position order as the foreach unit.

Not checked: suez_fit_models (LA mode), end-to-end package build, XL interleaved.
