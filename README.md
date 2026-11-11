# suez: an R package for mapping response-eQTLs

`suez` extends the `PANAMA` eQTL mapping framework (http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002330) to handle
- multiple conditions (e.g. stimulated/unstimulated), including response-eQTL mapping
- latent factor correction
- known kinship between individuals

## Installation 

Needs rstan / StanHeaders >= 2.26 and a C++17 compiler.

The easiest way to install is
```
if (!require("devtools")) install.packages("devtools", repos='http://cran.us.r-project.org')
devtools::install_github("ricardovialle/suez")
```

## Usage 

For example usage see [here](https://github.com/davidaknowles/dox/blob/523fb973e06daad34ce621df6a7644df27be02e2/code/panama_test.R). 

## Changes in this fork (v0.2)

Fixes found while validating suez against tensorQTL in simulations (repeated samples, relatedness,
within-individual G x condition). The models and function arguments are unchanged; output columns are added.

### Build and packaging

- Stan files use the `array[]` syntax (the old syntax was removed in Stan 2.33); `src/Makevars*` ask for
  C++17 (current StanHeaders). Needs rstan / StanHeaders >= 2.26.
- All dependencies are declared (`foreach`, `readr`, `stats`, `utils`; Bioconductor `VariantAnnotation`,
  `GenomicRanges`, `IRanges` as Suggests, only needed for a `TabixFile` genotype) and imported explicitly,
  which also removes the `magrittr::extract` / `rstan::extract` clash warning at load time.
- `R CMD check --no-manual` passes with one NOTE: `LICENSE` holds the Apache 2.0 text while DESCRIPTION
  says `GPL (>=3)` (inherited from upstream; to be resolved by the authors).

### `learn_covariance` (step 1)

- A numeric `condition` is regressed out as a slope. It was converted to a factor, so a continuous
  condition (one level per value) absorbed all the expression and every variance came out as 0.
- Variances and latent factors are fitted on the N - rank(X) contrasts orthogonal to the condition design
  (REML). The residuals of the regression have no variance along the design, and were fitted as
  N-dimensional data (an ML-type downward bias of the variances, ~ (N - k) / N).
- The `same_condition` component is dropped: after removing the condition effects it is identically zero
  (it was already excluded from `cov_suez`). `optimized_fit$par$s` therefore has no `same_condition` entry.
- Warns when `rstan::optimizing` does not report convergence.

### `map_interaction_qtl` (step 2)

- **Intercept-free likelihood**: the model is fitted on the N-1 contrasts orthogonal to the intercept
  (REML for the intercept only; LRTs of the other effects are unchanged). A GRM built from centred
  genotypes (plink, GCTA) has K 1 = 0: on the full data the intercept fits that direction exactly, the
  likelihood is unbounded as sigma2 -> 0, and the reported log-likelihoods were arbitrary (p-values off by
  up to ~4 orders of magnitude in a test with a 150-individual GRM).
- **Convergence**: Stan's `return_code` is checked (it was ignored); a failed fit gives NA instead of a
  silently wrong likelihood. Alternative models start from the default variances, because a null fit at the
  sigma2 = 0 boundary made the warm-started L-BFGS stop at its start (LR = 0), with a retry from the null
  variances if a fit lands below the null.
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

### Tests

`tests/` runs with `R CMD check` (about 15 s):

- `model_checks.R`: every `map_interaction_qtl` LRT against an exact profile-likelihood reference, with
  repeated samples, a centred GRM and a condition varying within individuals (including fits on the
  variance boundaries), effect sizes, bootstrap columns and the error paths.
- `learn_covariance_checks.R`: recovery of known noise, same-individual and latent-factor variances for a
  factor and a continuous condition.
