## Fitting helpers for map_interaction_qtl: the model y = X b + e,
## e ~ N(0, sigma2_k * K + sigma2 * I), fitted by maximum likelihood with
## stanmodels$suez_step_2 on data rotated by K's eigenvectors.
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
  stopifnot(nrow(K) == ncol(K), isSymmetric(unname(K), tol = 1e-8))
  n = nrow(K)
  Q = qr.Q(qr(matrix(1, n, 1)), complete = TRUE)[, -1, drop = FALSE]
  e = eigen(crossprod(Q, K %*% Q), symmetric = TRUE)
  if (min(e$values) < -1e-6 * max(abs(e$values)))
    stop("sample_kernel is not positive semi-definite (min eigenvalue ",
         signif(min(e$values), 3), ").", call. = FALSE)
  list(vectors = Q %*% e$vectors, values = pmax(e$values, 0))
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

# One ML fit; NULL unless Stan reports convergence.
fit_lmm = function(Uy, UX, lambda, init, debug = FALSE) {
  data = list(N = length(Uy), P = ncol(UX), U_transpose_x = UX, lambda = lambda, U_transpose_y = Uy)
  fit = tryCatch(stan_optimizing_wrapper(stanmodels$suez_step_2, data, init = init, as_vector = FALSE,
                                         verbose = debug),
                 error = function(e) { if (debug) message("Stan error: ", conditionMessage(e)); NULL })
  if (is.null(fit) || fit$return_code != 0 || !is.finite(fit$value)) return(NULL)
  fit
}

# Fit a model nested over `f0` (extra columns appended to UX). Starts from the
# default variances: a null fit at the sigma2 = 0 boundary stops L-BFGS at its
# starting point. Retries from the null fit's variances if that lands below the
# null likelihood; NULL if no fit reaches it (not a maximum).
fit_nested = function(Uy, UX, lambda, f0, debug = FALSE) {
  beta = as.array(c(f0$par$beta, numeric(ncol(UX) - length(f0$par$beta))))
  best = NULL
  for (init in list(list(sigma2 = 0.1, sigma2_k = 1.0, beta = beta),
                    list(sigma2 = f0$par$sigma2, sigma2_k = f0$par$sigma2_k, beta = beta))) {
    f = fit_lmm(Uy, UX, lambda, init, debug)
    if (!is.null(f) && (is.null(best) || f$value > best$value)) best = f
    if (!is.null(best) && best$value >= f0$value - 1e-6) break
  }
  if (is.null(best) || best$value < f0$value - 1e-3) return(NULL)
  best
}

# GLS estimate and SE of column k at the fitted variance components.
gls_coef = function(fit, UX, Uy, lambda, k = ncol(UX)) {
  o = fit$par$sigma2_k * lambda + fit$par$sigma2
  V = solve(crossprod(UX, UX / o))
  c(drop(V %*% crossprod(UX, Uy / o))[k], sqrt(V[k, k]))
}
