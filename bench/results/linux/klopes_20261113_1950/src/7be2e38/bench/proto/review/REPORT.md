# Independent reviewer — phase 1 report

(Written by the coordinator from the reviewer's hand-back; the reviewer could not write files outside its scripts.)

Files here: `validate.R` (compare two result tables, power/FDR vs truth, null calibration), `edge_cases.R`
(23-case suite for a candidate `fit_gene_proto`), `stress_exact.R` (global-optimum check on random kernels,
for methodological candidates), `interleave.R` (interleaved A/B timing), `memwatch.sh` (peak RSS summed over
the process tree), `baseline_checks.R`, `mutants/` (m1 one peak, m2 40-point grid, m3 LLT in grid), `out/`.

## Harness issues
| id | severity | issue | fix |
|---|---|---|---|
| H1 | High | `/usr/bin/time -l` maxrss = largest single process; doMC forks not summed (M_grm 4 cores: 307 MB reported vs 1076 MB tree) | `memwatch.sh`; compare RAM only at equal core counts |
| H2 | Medium | `proc.time()` misses CPU of forks not yet reaped | take user/sys from /usr/bin/time |
| H3 | High | timings noisy with concurrent agents (same run 9.4 vs 18.8 s); E vs P cores | interleaved prototype timing; serial authoritative runs on idle machine, >= 5 reps, median |
| H4 | High | size ladder uses default kernel with 2 reps -> eigenvalues {0, 2}: favours eigenvalue grouping | add L with GRM, unbalanced reps, a kernel learned by learn_covariance; report per kernel type |
| H5 | Medium | LD weak (max r2 0.48), no duplicate columns, fixed SNPs per gene | haplotype blocks, variable gene sizes |
| H6 | Medium | low signal: eGene power 0-23% | high-signal dataset for validation |
| H7 | Medium | failure paths (NA, df = 0) never exercised | edge_cases.R |
| H8 | Low | compare() rel beta divides by |beta| | divide by max(|beta|, se) |
| H9 | Low | compare sourceCpp prototypes with sourceCpp p0 | done in interleave.R |
| H10 | Low | TBB / OpenBLAS threads not capped | export RCPP_PARALLEL_NUM_THREADS |
| H11 | Info | ~220 MB fixed library cost in every RSS | also report RSS of a 0-gene run |

Verified: `prepare()` reproduces end-to-end inputs bitwise (bootstrap included); results bitwise identical across
reruns and 1/4/8 cores (S, M, M_dense, M_grm).

## Validation tool
```
Rscript bench/proto/review/validate.R compare <ref.rds> <cand.rds> [dataset]
Rscript bench/proto/review/validate.R truth   <res.rds> <dataset>
Rscript bench/proto/review/validate.R null    <lib> <dataset> [n_genes] [out.rds] [perm seed]
Rscript bench/proto/review/validate.R calib   <null_res.rds>
```
compare: per-column NA mismatches, max abs/rel diff, max |d -log10 p|, lead SNP and top-5 identity per gene,
per-gene Spearman, significant sets at 1e-3/1e-5/1e-8, eGenes (Bonferroni within gene, BH 5%), power/FDR, VERDICT.

Baseline: reruns and 4-core run identical (19,942 rows, EQUIVALENT). Null calibration (permuted genotypes):
S p_geno lambda_GC 1.03 (KS p 0.23), p_interact 1.01; M_grm 0.98/0.98; bootstrap LR slightly heavier than
chi2(1) on S (finite-sample ML LRT, a property of the baseline).

## Numerical fragility (what candidates must preserve)
1. Grid uses normal equations; final values/derivatives/beta use direct residuals (Profile). Moving the final
   values to normal equations loses accuracy (weights up to 1e5 on the grid, 1e12 in Newton).
2. +-Inf sentinels (h = 1 with zero eigenvalues -> -Inf; non-positive variance -> +Inf -> NA). NaN would change
   peak detection and NA patterns.
3. Rank-deficient designs: LDLT pseudo-solves; beta/se via LU inverse (NaN for collinear covariates — a
   pre-existing issue). LLT in the grid changed l by 0.33.
4. Mixed df within a gene (absent G x C column).
5. Grid/peaks are methodology: MAX_PEAKS = 3, 80 points on [-5, 5] + h = 0, 1; Newton |x| <= 12.
6. Baseline within 1e-8 of an exact optimiser: the 1e-6 tolerance is meaningful.
7. ok_nested threshold 1e-6.
8. Bootstrap depends on eigenvector signs from eigen(): any change in how U is computed changes the draws
   (judge distributionally); threaded candidates must draw normals sequentially first.
9. Empty paths: p0 = 0, q = 0, one SNP, small odd n.

Key caveat: passing tolerances does not prove a pure speed-up — the 40-point-grid mutant is a methodological
change yet passes everything (0.69x time). Methodological status must come from reading the diff.

Checklist per candidate: read diff; edge_cases.R; validate.R compare on S, M_grm, M_learn, M_fac3, M_eqtl,
M_cov (+ LA if grid/refine change); null + calib; interleave.R timing; memwatch.sh RAM at equal cores.
