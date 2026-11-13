# Global-optimum stress test for methodological candidates (grid density, peaks refined, tolerances):
# random 3-block kernel spectra (often bimodal likelihoods in h, narrow peaks, boundary optima); l0 of the
# baseline and of a candidate vs the exact ML (1601-point log-delta grid, every local max refined).
#   Rscript bench/proto/review/stress_exact.R <candidate.cpp> [n_cases]
source("bench/proto/harness.R")
a = commandArgs(TRUE); Rcpp::sourceCpp(a[1]); nc = if (length(a) > 1) as.integer(a[2]) else 300
exact = function(UX, Uy, lam) {
  f = function(x) { h = 1 / (1 + 10^-x); d = h * lam + (1 - h); W = 1 / d
    b = solve(crossprod(UX, UX * W), crossprod(UX, Uy * W))
    -0.5 * length(Uy) - 0.5 * sum(log(sum((Uy - UX %*% b)^2 * W) / length(Uy) * d)) }
  xs = seq(-8, 8, length.out = 1601); v = sapply(xs, f); pk = which(diff(sign(diff(v))) == -2) + 1
  c(max(c(v[c(1, 1601)], vapply(pk, function(i) optimize(f, xs[c(i - 1, i + 1)], maximum = TRUE, tol = 1e-10)$objective, 0))), length(pk))
}
set.seed(99); n = 120
r = t(replicate(nc, {
  Q = qr.Q(qr(matrix(1, n, 1)), complete = TRUE)[, -1]
  m = sort(sample(10:90, 2)); sz = c(m[1], m[2] - m[1], n - 1 - m[2])
  lam = rep(10^runif(3, -4, 3), sz) * exp(rnorm(n - 1, 0, 0.3))        # 3 eigenvalue blocks, random scales
  z = rnorm(n - 1) * sqrt(rep(10^runif(3, -2, 2), sz))                 # random variance per block
  K = Q %*% (lam * t(Q)) + 1 / n; e = suez:::kernel_contrast_eigen(K)
  y = drop(Q %*% z); cond = rnorm(n); UX0 = crossprod(e$vectors, (cond - mean(cond)) / sd(cond))
  Uy = drop(crossprod(e$vectors, y)); G = cbind(rbinom(n, 2, 0.3), rbinom(n, 2, 0.4))
  ex = exact(UX0, Uy, e$values)
  b = baseline_fit(e$vectors, e$values, Uy, UX0, G, matrix(0, n, 0), FALSE)$l0
  c_ = fit_gene_proto(e$vectors, e$values, Uy, UX0, G, matrix(0, n, 0), FALSE)$l0
  c(peaks = ex[2], base_minus_exact = b - ex[1], cand_minus_exact = c_ - ex[1], cand_minus_base = c_ - b)
}))
cat(sprintf("%d cases (%d with >= 2 interior peaks)\n", nc, sum(r[, 1] >= 2)))
for (j in 2:4) cat(sprintf("%-17s below by > 1e-6: %3d   max shortfall %.3g   max excess %.3g\n", colnames(r)[j],
                           sum(r[, j] < -1e-6), -min(r[, j]), max(r[, j])))
