# Round 12 — independent reviewer / benchmark specialist (transcribed by the coordinator from the hand-back)

Files: validate12.R (one-command suite), run12.R (runner incl. permutation / null runs), edge12.R, stress12.R,
datasets12.R + data/ (P_many10, XL_grm_many, R_real), flags_check.R/.log, mem_gene.R, results/ (verdict TSVs),
out/ (cached runs), pkg_mut_1peak/ (negative control, bench/lib/review_mut_1peak, MAX_PEAKS = 1).

## 1. Audit
| id | sev | finding | fix |
|---|---|---|---|
| A1 | High | P_many is not 3000 x 10: spacing 3e4 with cisdist 1e5 gives ~69 SNPs/gene, 206,669 tests | P_many10 (spacing 2.5e5, 10 SNPs/gene, 29,986 tests) |
| A2 | High | no large-N workload with many genes (XXL wall ~75% set-up) | XL_grm_many (N = 2000, GRM, 400 genes x 25 SNPs); report set-up vs per-gene separately |
| A3 | High | equal SNPs/gene everywhere; GRM never with covariates / missing / unbalanced — the slow path (~2.4 ms/test vs 0.12 on M_cov) | R_real: N = 850 (400 ind., 1-3 samples, condition within individual), GRM, 10 covariates, 1% missing, lognormal 0-1500 SNPs/region (8% empty), clustered genes, rare SNPs |
| A4 | Med | wide windows / >1M-row tables never benchmarked | mem_gene.R; propose P_big (20k genes x 100 SNPs) |
| A5 | Med | fixed implementation order per repetition | alternate order |
| A6 | Low | cpu_util mixes run.R wall with /usr/bin/time user+sys (includes start-up) | use run.R user/sys; report a 0-gene RSS |
| A7 | Low | bench.sh set -e skips diagnostics on failure; temp logs leak | grep ... \|\| true |
| A8 | Low | MKL/FLEXIBLAS/GOTO threads not pinned | pin in bench.sh |
| A9 | Med | harness.R defaults to bench/lib/baseline (0.4) | default base12 |
| A10 | Low | prepare(): fixed cis window, no LA; compare() beta / abs(beta) | use validate12.R |
| A11 | Info | sourceCpp lacks -DEIGEN_NO_DEBUG; measured ratio vs installed 0.93-1.01, rep spread +-8% | claims < 10% need >= 5 reps |
| A12 | Med | round-1 tools: validate.R Inf/Inf pair skips a column's tolerance check; edge_cases.R passes Inf vs finite; LA columns unchecked; lead-SNP ties | fixed in validate12.R / edge12.R |
| A13 | Info | profile12.R appends duplicate rows on rerun | - |
Not covered: maf_threshold > 0, checkpoint_dir, TabixFile input, NA expression, cisdist = 1e6.

## 2. Suite
`Rscript bench/proto/r12/review/validate12.R <cand> [--ref base12] [--ds all|quick|list] [--parts equiv,stat,repro,edge,stress] [--boot exact|dist] [--jobs 2] [--rerun]`
-> results/v12_<cand>_vs_<ref>.tsv (cand, ref, part, dataset, check, value, threshold, status, detail), VERDICT line.
equiv: rows/cols/NA/Inf pattern; l 1e-9; boot l 1e-5 (--boot exact); beta (rel. max(|beta|, se)) and se 1e-6 incl.
LA columns; df identical; other numeric 1e-12; per p column: lead SNP, significant sets 1e-3/1e-5/1e-8, eGenes
(Bonferroni within gene then Bonferroni / BH 5%), within-gene rank agreement (1e-6 relative band; borderline
reported, not failed); bootstrap LR and l_boot_geno - l_geno KS + quantiles (gate only with --boot dist).
stat: power/FDR vs truth; null lambda_GC difference vs reference (FAIL > 0.05, WARN > 0.02).
repro: identical() on S, M_dense, M_grm, LA for a rerun and 2 / 4 cores. edge: 23 edge_cases inputs (strict
NA/Inf; l 1e-9 PASS, 1e-6 WARN). stress: 150 random 3-block spectra + bimodal cases vs the exact optimum.

## 3. Sanity
base12 vs itself: 1105 checks PASS (15 min). final vs base12: 1105 PASS (diffs only L_grm, XL_grm, R_real from
round-11 rotation: l <= 3.2e-12, beta/se <= 4.1e-8, boot <= 2.5e-6). Negative control MAX_PEAKS = 1: FAIL on
M_learn (l_interact 0.0067), edge bimodal (0.067), stress bimodal (1.05e-5); 11 of 12 quick datasets PASS =>
edge and stress parts are mandatory for every candidate.

## 4. Risks
R1 (High, memory) per-gene memory O(N S) unblocked (G, Rcpp copy, Gs, Is[j], RG): N = 2000, 4000 SNPs -> 883 MB
peak (7.4x the genotype matrix); fix: lazy rotation in 256-column chunks aligned with SNP blocks.
R2 (Med) fit_models maps cols[j] through a temporary NumericMatrix: use-after-free if an integer matrix is passed.
R3 (Med) failing genes dropped silently (.errorhandling = 'remove').
R4 (Med, perf) GRM + covariates: Profile/Newton costs R x (p0+k+1)^2/2 pairs per evaluation without the grid's
Schur shortcut; R_real ~20x M_cov per test.
R5 (Low) U'y per gene in R is O(N^2) (one pass over U per gene).
R6 (Low) per-gene constant set-up (Spectrum, GridW, x0 moments/grid) cacheable (~55 us/gene).
R7 (Low) rotate() rebuilds UZ per chunk/term. R8 kernel_contrast_eigen holds 5-6 N x N at peak.
R9 exactness stress-tested only at n <= 400. R10 bootstrap draws depend on the eigen basis (Linux --boot dist).

## 5. Kernel-scale / boundary gap (p_prof2) and the `tests` part
Gap: p_prof2 started Newton from -Inf; Newton's x = log10(h/(1-h)) is capped at +-12, so a grid optimum at h = 0 / 1
could not be reached (indistinguishable at f = 1; up to 1e-5 below the ML at f = 1e6). Only tests/model_checks.R
caught it: no suite case combined a boundary optimum with a large kernel scale.
kscale12.R (inside --parts edge): case kernel_scale (model_checks.R data, sample_kernel = f K, f = 1, 1e6, 1e-6; 144
fits per scale vs an exact optimiser over log10 delta in [-25, 25] + delta = 0 / Inf); gates f = 1e6 within 1e-9 of
exact, invariance within 1e-6, f = 1 / 1e-6 within 1e-9 or no worse than the reference + 1e-9 (base12 is 9.0e-8 below
exact at f = 1e-6, pre-existing). Cases boundary_h0 / boundary_h1 (exact ML at h = 0 / h = 1; scales 1, 1e3, 1e6)
gate 1e-9. `tests` part: every tests/*.R of the candidate's source against the candidate lib (FAIL on FAIL line,
non-zero exit, or no PASS line); source from --src or bench/proto/r12/<role>/pkg_<name>; now in the default parts.
Verification (--parts edge,tests): base12 PASS 36/36; integ_c3 FAIL (5): kernel_scale f=1e6 1.16e-5, boundary_h0
f=1e3 3.9e-8 / f=1e6 3.9e-5, invariance 1.16e-5, model_checks FAIL; numerical_comb FAIL (same); integ_c3b PASS.
Full run integ_c3b vs base12 (--ds all --boot dist --jobs 2): PASS with warnings, 1118 checks (977 PASS, 3 WARN, 0
FAIL, 138 INFO), 334 s. WARN: L bootstrap KS p 0.0096 (smallest of 29 KS p; P(min <= that) ~ 0.24; quantiles agree);
edge collinear_interaction_and_snps 2 NA/Inf mismatches (report-only; same on integ_c3 / numerical_comb); edge
near_perfect_fit 7.95e-7 (report-only). Max l 1.15e-10 (XXL), beta 1.2e-7, se 3.0e-8 relative; leads, significant
sets, eGenes, ranks 0 failures; stat 34/34; repro identical 1/2/4 cores; stress 0 misses; tests 6/12/17 PASS.
Cost +4 min per candidate. Not covered: boundary optima at large N, kernel scale with covariates, scales > 1e6.
