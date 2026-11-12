# Changelog

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
  h = sigma2_k / (sigma2_k + sigma2), with the effects and total variance in closed form (grid plus
  golden-section search, boundaries included). This replaces `rstan::optimizing`, which ignored failed
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
- `maf_threshold` filters on the minor allele frequency. It used the ALT allele frequency, so SNPs with a
  common ALT allele (e.g. ALT frequency 0.98) were kept whatever the threshold.
- `checkpoint_dir` files are written inside the directory (`dir/gene.txt.gz`; without a trailing slash they
  went next to it as `dirgene.txt.gz`).

#### Performance

- Genes run in parallel with the registered `foreach` backend (e.g. `doMC::registerDoMC(8)`); SNPs within a
  gene run in sequence. Each gene gets its own seed drawn from the caller's RNG, so with `set.seed()` the
  bootstrap draws are the same for any number of cores.
- Step 2 is ~3x faster on one core (C++ fits); 64 genes x 20 SNPs with bootstrap took 18 s on 1 core and
  1.2 s on 8 (Apple M5, 4 performance + 6 efficiency cores).
- Step 1 time no longer grows with the number of genes (see above).
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
