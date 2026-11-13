# Numerical computing specialist: report
See the final hand-back message (same content). Summary:

- **Recommended: p7_exact.cpp.** It keeps the method unchanged (grid, peaks, Newton, RNG order) and changes how
  the sums are computed:
  - GEMM grid with a Schur complement on the X0 block;
  - eigenvalue-group moments;
  - rotation that skips duplicate sample rows;
  - plain-loop small solves;
  - explicit-residual fallback in Profile;
  - LU beta/se, as in p0.
- **Speed-up of p7 over sourceCpp p0** (median of 3 interleaved reps): S 8.8x, S_rare 7.9x, M_fac3 13x,
  M_cov 8.3x, L 9.9x, M_unbal 10.8x, M_learn 8.8x, M_lc 9.1x, M_grm 4.2x, M_num 4.3x, L_grm 4.2x, M_eqtl 2.4x.
- **Correctness:**
  - all log-likelihoods are within 1e-6 (bootstrap columns ≤ 4.5e-7, others ≤ 1e-12);
  - 0 NA mismatches and identical df on 12 datasets;
  - |Δbeta| ≤ 5e-8 se;
  - edge_cases.R: 21 PASS, 0 FAIL, 2 DIFF on report-only cases (duplicated_covariate, near_perfect_fit).
- **Rejected:** p6 (single-product Profile, mixed ±), p8 (frexp logdet, 0.67x), -O3 / -mcpu=native (≤ 3%).
- Raw data: p7_all.txt. Scripts: run.R, tint.R, betadiag.R, gridtest.R, edge_detail.R, timed.R.
