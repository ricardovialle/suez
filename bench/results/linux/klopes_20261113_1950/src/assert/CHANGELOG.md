# Changelog

## Unreleased (round 11: BLAS/LAPACK and large-N review)

- `kernel_contrast_eigen` picks its eigen solver from the BLAS/LAPACK R links to (`options(suez.eigen_solver =
  "auto" | "eigen" | "lapack")`): 0.5 always used Eigen's solver, measured only against R's reference LAPACK
  (where it is 2x faster). With OpenBLAS, LAPACK is 4.4x faster at N = 2000 and 6.6x at N = 5000 (one thread;
  38 -> 6 s). Data fits identical; **bootstrap draws for a given seed change** where LAPACK is now used.
- `suez_fit_gene` / `suez_fit_models` rotate all SNPs of a gene at once (256-column chunks) instead of per
  moment block: at large N the blocks are 16 SNPs wide and each re-read the N x N-1 eigenvector matrix (200 MB at
  N = 5000), which made 0.5 slower than 0.3 there. Same results to rounding (<= 3e-11).
- Non-finite genotype or term values are refused in C++ (they broke `std::sort`'s ordering: undefined
  behaviour, a possible crash with GCC's libstdc++). `map_interaction_qtl` already imputes; direct callers of
  `suez_fit_gene` (e.g. eqtl_smk) are protected now.
- `bench/linux/run_checks.sh`: build + tests (also with libstdc++ assertions) + eigen and fit timings vs 0.3/0.5.

## 0.5

Performance: step 2 is 1.9-7.4x faster on one core and uses 22-41% less peak RAM; step 1 is 1.6-2.0x faster
with 11-32% less RAM. Same models, grid, peak refinement and Newton steps; log-likelihoods, p-values,
effect sizes and SEs agree with 0.4 to 1e-9 or better (18 benchmark datasets, `bench/`). Two documented
differences: the parametric bootstrap draws for a given seed differ (same distribution), and step 1 converges
along a slightly different path (~1e-4 relative in `cov_suez`).

#### `map_interaction_qtl` (step 2)

- The profile likelihood is evaluated from moments of contrasts grouped by kernel eigenvalue: tied non-null
  eigenvalues (within 1e-12 relative) share a weight; the null eigenvalues are kept exact by a Taylor
  expansion in their rounding-level deviation, checked at every h, with a per-contrast fallback. The grid over
  h is one matrix product for all grid points and SNPs; Newton uses the same moments and switches to direct
  residuals where they would cancel. SNPs are processed in blocks of ~4 MB of moments (memory does not grow
  with covariates x SNPs), and repeated sample rows are rotated once. Local-ancestry models use the same
  engine (~5x faster).
- `kernel_contrast_eigen`: the contrast basis is applied as a Householder reflection (O(N^2) instead of dense
  O(N^3) products) and the eigen-decomposition uses Eigen's symmetric solver: 23 -> 2.5 s at N = 2000. Inside
  tied eigenvalues the basis differs from `eigen()`'s, so **bootstrap draws (`l_boot_*`) for a given
  `set.seed()` differ from 0.4**; their distribution is the same (KS tests on 17 datasets) and results remain
  identical across numbers of cores.
- Wald estimates (`beta_*`, `se_*`) are computed on the identified columns: a design column that is a
  combination of the others (1 - R^2 < 1e-14, e.g. local-ancestry counts when an ancestry is absent at the
  locus) is dropped and its estimate is NA. 0.4 inverted the singular matrix (NaN SEs for identified terms).
- `easy_impute` reconstructs missing genotypes from the eigen-decomposition of the smaller Gram matrix (2.2x
  faster, identical imputed genotypes).
- `U` is passed to C++ without a copy per gene; the default kernel is built as Z Z' and freed after the
  eigen-decomposition; `quantile_normalize` works row by row (input + 1 copy instead of ~6).

#### `learn_covariance` (step 1)

- The Stan model sums the fixed components as one matrix-vector product and the latent factors as one
  `tcrossprod`: fewer autodiff nodes, 32-60% less peak RAM at N = 400-1000 and 1.6-2.0x faster. Same
  objective; the optimiser path differs (`cov_suez` ~1e-4 relative, step-2 -log10 p within 0.003).

#### Package

- rstan is no longer attached at load time (`learn_covariance` calls `rstan::optimizing`): `library(suez)`
  takes ~130 MB less RAM and ~0.6 s less, and each forked worker ~70 MB less.

#### Tests and benchmarks

- New tests: likelihoods invariant to the kernel scale (x1e6, x1e-6), and an absent third ancestry gives the
  two-ancestry estimates.
- `bench/`: simulated workloads, benchmark harness, specialist prototypes, experiment log and results
  (`bench/REPORT.md`).

## 0.4

Local-ancestry (Tractor-style) eQTL and response-eQTL mapping for admixed samples.

#### `map_interaction_qtl` (step 2)

- **`local_ancestry`**: per-ancestry ALT dosages and haplotype counts (Tractor's `extract_tracts.py`
  output, read with the new `read_tractor()`) replace `genotype`. Each SNP is fitted with the models of
  Tractor (Atkinson et al. 2021) inside the suez mixed model: local ancestry plus one genotype term per
  ancestry, and their interactions with the condition. New tests: any ancestry-specific eQTL (`p_geno`, k
  df), eQTL effects differ between ancestries (`p_het`), any ancestry-specific response (`p_interact`, with
  local ancestry x condition in the null), responses differ (`p_interact_het`), and per-ancestry Wald
  estimates (`beta_<anc>`, `p_<anc>`, `beta_interact_<anc>`, ...). Works for any number of ancestries;
  an ancestry with fewer than 5 copies of either allele at a SNP is dropped for that SNP. Bootstrap under the
  no-response model.
- With one sample per individual and the default kernel, the per-ancestry betas, SEs and p-values equal
  Tractor's `run_tractor.R` linear model (Wald SEs on the residual degrees of freedom with a t reference),
  also on the Tractor tutorial data (61 ASW individuals, chr22, four linear phenotypes).
- **`covariates`**: fixed covariates (e.g. global ancestry, PCs, age, sex) added to every model, in both
  modes.
- A `condition` with a single value maps eQTLs without interaction tests (it was an error for a factor).
- `src/lmm.cpp`: a general fitter for any set of per-SNP terms (`suez_fit_models`), using the same
  batched grid and Newton search.

#### Tests and documentation

- `tests/local_ancestry_checks.R`: every local-ancestry likelihood against the exact profile-likelihood
  reference, equality with Tractor's linear model, three ancestries, dropped rare ancestries and
  `read_tractor()`.
- README: local-ancestry section (inputs, models, output columns, recommendations) and an updated
  schematic.
- Vignette section 8: local-ancestry mapping on a simulated admixed cohort (`inst/extdata/admixed`, Tractor
  file format, generated by `data-raw/make_example_data.R`).
- Simulations behind the recommendations (300 admixed individuals, 2 conditions, sibling pairs, LD that
  differs by ancestry): a GRM in step 1 does not absorb ancestry-specific effects (same estimates and
  power) and removes the inflation from relatives; a response present only in the minor ancestry is found
  far more often than by the pooled test; LD that differs by ancestry makes tag SNPs look ancestry-specific.

## 0.3

Faster and more accurate search over h in step 2; results equal 0.2 except where 0.2 missed the maximum.

#### `map_interaction_qtl` (step 2)

- `maf_threshold` filters on the minor allele frequency. It used the ALT allele frequency, so SNPs with a
  common ALT allele (e.g. ALT frequency 0.98) were kept whatever the threshold.

#### Performance

- Step 2 fits all cis-SNPs of a gene in one C++ call: the h grid is evaluated for every SNP together
  (condition terms once per grid point, genotype terms as matrix products across SNPs), then each SNP is
  refined with Brent's method; the grid is dense near h = 1 and the refinement there works in
  log10(1 - h), where a kernel with null directions can put the maximum within 1e-6 of the boundary.
  Results match the previous version (log-likelihoods within 1e-13); 4.7-5.6x faster per core
  (e.g. 731 -> 4,120 SNPs/s on one core, 400 genes x 50 SNPs, 200 samples, bootstrap).
- The search over h follows GEMMA / FaST-LMM: a grid uniform in log10(delta), delta = sigma2_k / sigma2
  (80 points on [-5, 5] plus h = 0 and h = 1), then safeguarded Newton steps with analytic first and second
  derivatives around every local maximum of the grid, keeping the best. A grid uniform in h missed narrow
  or second peaks when the kernel's eigenvalues span orders of magnitude (e.g. a learned `cov_suez`): on the
  vignette data 4 of 13,801 fits were up to 0.065 log-likelihood below the maximum. Checked against a
  1,000-point reference on seven datasets: no fit below it.
- cis-SNPs are found with a per-chromosome index of `snploc` sorted by position, built once (it was a scan
  of all of `snploc` per gene, ~20 ms per gene per million SNPs). SNPs in `snploc` missing from `genotype`
  are skipped (they made the gene fail).
- Genes without missing genotypes skip the SVD imputation, and results are combined once at the end
  instead of pairwise. With 1 million SNPs in `snploc`: 2,957 -> 4,719 tests/s on one core, 11,766 ->
  18,368 on 8 cores; results equal the previous version except where it missed the maximum.

#### Tests and documentation

- `model_checks.R`: a likelihood with two near-equal peaks in h (fails if only the best grid point is
  refined), and `maf_threshold` with the ALT allele coded either way.
- The vignette sets a seed before `learn_covariance` (rstan draws its seed from R's RNG), so it renders
  the same every time.

## 0.2 (fork: ricardovialle/suez)

Fixes found while validating suez against tensorQTL in simulations (repeated samples, relatedness,
within-individual G x condition). The models and function arguments are unchanged; output columns are added.

#### Build and packaging

- Stan files use the `array[]` syntax (the old syntax was removed in Stan 2.33); `src/Makevars*` ask for
  C++17 (current StanHeaders). Needs rstan / StanHeaders >= 2.26.
- All dependencies are declared (`foreach`, `readr`, `stats`, `utils`; Bioconductor `VariantAnnotation`,
  `GenomicRanges`, `IRanges` as Suggests, only needed for a `TabixFile` genotype) and imported explicitly,
  which also removes the `magrittr::extract` / `rstan::extract` clash warning at load time.
- License: GPL (>= 3), as in DESCRIPTION; `LICENSE` now holds the GPL v3 text (it held Apache 2.0, which
  conflicted with DESCRIPTION and with the GPL v3 code from rstanarm/rstantools in the package).
- `R CMD check --no-manual`: Status OK.

#### `learn_covariance` (step 1)

- A numeric `condition` is regressed out as a slope. It was converted to a factor, so a continuous
  condition (one level per value) absorbed all the expression and every variance came out as 0.
- Variances and latent factors are fitted on the N - rank(X) contrasts orthogonal to the condition design
  (REML). The residuals of the regression have no variance along the design, and were fitted as
  N-dimensional data (an ML-type downward bias of the variances, ~ (N - k) / N).
- The `same_condition` component is dropped: after removing the condition effects it is identically zero
  (it was already excluded from `cov_suez`). `optimized_fit$par$s` therefore has no `same_condition` entry.
- Warns when `rstan::optimizing` does not report convergence.
- The Stan model takes the sufficient statistic `sum_g y_g y_g'` (samples x samples) instead of looping over
  genes, so run time no longer grows with the number of genes (same likelihood; ~10-50x faster in tests).

#### `map_interaction_qtl` (step 2)

- **Intercept-free likelihood**: the model is fitted on the N-1 contrasts orthogonal to the intercept
  (REML for the intercept only; LRTs of the other effects are unchanged). A GRM built from centred
  genotypes (plink, GCTA) has K 1 = 0: on the full data the intercept fits that direction exactly, the
  likelihood is unbounded as sigma2 -> 0, and the reported log-likelihoods were arbitrary (p-values off by
  up to ~4 orders of magnitude in a test with a 150-individual GRM).
- **Exact fits**: each model is fitted in C++ (`src/lmm.cpp`) by maximising the likelihood over
  h = sigma2_k / (sigma2_k + sigma2), with the effects and total variance in closed form (grid plus a
  one-dimensional search, boundaries included; the search was changed in 0.3). This replaces `rstan::optimizing`, which ignored failed
  fits, and whose L-BFGS stopped short near sigma2 = 0 (LR = 0 when warm-started on the boundary; ~0.03
  log-likelihood short otherwise). A fit that is not finite gives NA.
- **Scaling**: covariate, genotype and interaction columns are standardized before fitting (invariant
  likelihood; L-BFGS stopped short on raw-scale conditions such as age).
- **Checks**: stops on a kernel that is not positive semi-definite, or when the condition design fits every
  null direction of the kernel (e.g. every within-individual contrast with the default kernel), where sigma2
  is not identifiable.
- **Outputs**: `p_geno`, `p_interact` with `df_interact` (1 for a numeric condition, levels - 1 for a
  factor), and the GLS genotype effect `beta_geno` / `se_geno` per ALT allele.
- **Bootstrap**: drawn in the rotated space from the fitted variances (no Cholesky of a possibly singular
  covariance), with a null fit on the bootstrap data.
- The `/dev/null` connection opened when `debug = FALSE` is closed on exit (it leaked one per call).
- `checkpoint_dir` files are written inside the directory (`dir/gene.txt.gz`; without a trailing slash they
  went next to it as `dirgene.txt.gz`).

#### Performance

- Genes run in parallel with the registered `foreach` backend (e.g. `doMC::registerDoMC(8)`); SNPs within a
  gene run in sequence. Each gene gets its own seed drawn from the caller's RNG, so with `set.seed()` the
  bootstrap draws are the same for any number of cores.
- Step 2 is ~3x faster on one core (C++ fits); 64 genes x 20 SNPs with bootstrap took 18 s on 1 core and
  1.2 s on 8 (Apple M5, 4 performance + 6 efficiency cores).
- Step 1 time no longer grows with the number of genes (see above).

#### Tests

`tests/` run with `R CMD check` (about 10 s):

- `model_checks.R`: every `map_interaction_qtl` LRT against an exact profile-likelihood reference, with
  repeated samples, a centred GRM and a condition varying within individuals (including fits on the
  variance boundaries), effect sizes, bootstrap columns, identical results on 1 and 2 cores,
  checkpoint files and the error paths.
- `learn_covariance_checks.R`: recovery of known noise, same-individual and latent-factor variances for a
  factor and a continuous condition.

#### Documentation

- README rewritten: installation with `pak` (`devtools::install_github()` is deprecated), a runnable example,
  input formats, output columns, checkpoint files, bootstrap calibration and parallel use. Change notes moved
  here.
- Vignette `vignettes/suez.Rmd` (rendered at https://ricardovialle.github.io/suez/doc/suez.html): raw counts and a VCF to eQTL calls, with
  simulated example data in `inst/extdata` (generated by `data-raw/make_example_data.R`).

## 0.1

- Original release by David A. Knowles (2017).
