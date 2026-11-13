## Fitting helpers for map_interaction_qtl: the model y = X b + e,
## e ~ N(0, sigma2_k * K + sigma2 * I), fitted by maximum likelihood
## (src/lmm.cpp, all SNPs of a gene per call) on data rotated by K's eigenvectors.
##
## The likelihood is that of the N-1 contrasts orthogonal to the intercept
## (REML for the intercept only). Every model has an intercept, so likelihood-
## ratio tests of the other effects are unchanged. This matters for a kernel
## with K 1 = 0, e.g. a GRM built from centred genotypes (plink, GCTA): on the
## full data the intercept fits that direction exactly and the likelihood grows
## without bound as sigma2 -> 0, so the reported log-likelihoods are arbitrary.

# Eigen-decomposition of K restricted to the contrasts: vectors = Q V
# (N x N-1, Q an orthonormal basis of the complement of 1). Refuses a kernel
# that is not positive semi-definite.
kernel_contrast_eigen = function(K) {
  K = as.matrix(K)
  stopifnot(nrow(K) == ncol(K), is_symmetric(K, tol = 1e-8))
  n = nrow(K)
  # = qr.Q(qr(1), complete = TRUE)[, -1] (same Householder product, column by column), without the N x N
  # identity and the complete Q
  I = matrix(0, n, n - 1); I[cbind(seq_len(n - 1) + 1, seq_len(n - 1))] = 1
  Q = qr.qy(qr(matrix(1, n, 1)), I); rm(I)
  e = eigen(crossprod(Q, K %*% Q), symmetric = TRUE)
  if (min(e$values) < -1e-6 * max(abs(e$values)))
    stop("sample_kernel is not positive semi-definite (min eigenvalue ",
         signif(min(e$values), 3), ").", call. = FALSE)
  list(vectors = Q %*% e$vectors, values = pmax(e$values, 0))
}

# isSymmetric(unname(K), tol) (all.equal: mean relative difference over the unequal entries), one column
# at a time: isSymmetric holds ~4 N x N temporaries.
is_symmetric = function(K, tol) {
  d = a = m = 0
  for (j in seq_len(ncol(K))) {
    x = K[, j]; y = K[j, ]; u = x != y
    if (anyNA(u)) return(FALSE)
    d = d + sum(abs(x[u] - y[u])); a = a + sum(abs(x[u])); m = m + sum(u)
  }
  m == 0 || (if (a / m > tol) d / a else d / m) <= tol
}

# Null directions of the kernel (eigenvalue ~ 0, e.g. contrasts between samples
# of the same individual with the default kernel) only have variance sigma2. If
# the fixed effects fit all of them exactly the likelihood is unbounded as
# sigma2 -> 0. `need` = columns the alternative models add that can vary within
# those directions.
check_identifiable = function(eig, X0, need = 1) {
  null = eig$values <= 1e-8 * max(eig$values)
  m = sum(null)
  if (m == 0) return(invisible(Inf))
  d = if (ncol(X0)) svd(crossprod(eig$vectors[, null, drop = FALSE], X0), 0, 0)$d else 0
  r = sum(d > 1e-8 * max(d, 1))   # svd rank: qr's rank would count a ~0 column
  if (m - r < need) stop(
    "the sample kernel leaves ", m, " direction(s) with no kernel variance (e.g. contrasts ",
    "between samples of the same individual); the condition design fits ", r, " and the ",
    "tested effects need ", need, " more, so sigma2 is not identifiable.", call. = FALSE)
  invisible(m - r)
}

# Centre and scale columns (the likelihood is invariant; L-BFGS stops short on
# raw-scale columns such as age).
standardize_cols = function(M) {
  if (ncol(M)) {
    M = scale(M)
    attr(M, "scaled:center") = attr(M, "scaled:scale") = NULL
  }
  M
}
