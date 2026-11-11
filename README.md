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

Fixes to `map_interaction_qtl`, found while validating suez against tensorQTL in simulations
(repeated samples, relatedness, within-individual G x condition). The model and the API are unchanged;
output columns are added.

- **Build**: Stan files use the `array[]` syntax (the old syntax was removed in Stan 2.33) and
  `src/Makevars*` ask for C++17 (current StanHeaders).
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

`tests/model_checks.R` checks every LRT against an exact profile-likelihood reference, including fits on
the variance boundaries, plus the error paths. Not changed: `learn_covariance` (step 1).
