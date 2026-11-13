# Bimodal-weighted stress test: random 3-block spectra tuned (uniroot on the third block's scale) so that
# the likelihood in h has two peaks whose heights differ by a random gap (+-1e-5 .. 1e-1); l0 of the baseline
# and each candidate vs the exact ML.   Rscript bench/proto/review/stress_bimodal.R [n_cases] cand.cpp...
source("bench/proto/harness.R")
a = commandArgs(TRUE); nc = as.integer(a[1]); impl = list(baseline = baseline_fit)
for (f in a[-1]) { e = new.env(); Rcpp::sourceCpp(f, env = e); impl[[sub("\\.cpp$", "", basename(f))]] = e$fit_gene_proto }
peaks = function(UX, Uy, lam) {
  f = function(x) { h = 1 / (1 + 10^-x); d = h * lam + (1 - h); W = 1 / d
    b = solve(crossprod(UX, UX * W), crossprod(UX, Uy * W)); -0.5 * length(Uy) - 0.5 * sum(log(sum((Uy - UX %*% b)^2 * W) / length(Uy) * d)) }
  xs = seq(-8, 8, length.out = 801); v = sapply(xs, f); pk = which(diff(sign(diff(v))) == -2) + 1
  vapply(pk, function(i) optimize(f, xs[c(i - 1, i + 1)], maximum = TRUE, tol = 1e-11)$objective, 0)
}
set.seed(17); n = 120; Q = qr.Q(qr(matrix(1, n, 1)), complete = TRUE)[, -1]; sz = c(40, 40, 39)
rows = list(); tries = 0
while (length(rows) < nc && tries < 20 * nc) {
  tries = tries + 1
  lam = rep(10^c(runif(1, -3.5, -2.5), runif(1, -0.5, 0.5), runif(1, 1, 2)), sz) * exp(rnorm(n - 1, 0, 0.3))
  K = Q %*% (lam * t(Q)) + 1 / n; e = suez:::kernel_contrast_eigen(K); cond = rnorm(n)
  UX0 = crossprod(e$vectors, (cond - mean(cond)) / sd(cond))
  z = rnorm(n - 1) * sqrt(rep(c(runif(2, 0.01, 2), runif(1, 0.01, 60)), sz))
  y_of = function(t) drop(Q %*% (z * rep(c(1, 1, t), sz)))
  gap = function(t) { p = peaks(UX0, drop(crossprod(e$vectors, y_of(t))), e$values); if (length(p) < 2) NA else p[1] - p[length(p)] }
  ts = 10^seq(-0.5, 0.5, length.out = 9); g = sapply(ts, gap)
  j = which(!is.na(g[-1]) & !is.na(g[-9]) & sign(g[-1]) != sign(g[-9]))[1]; if (is.na(j)) next
  target = sample(c(-1, 1), 1) * 10^runif(1, -5, -2)
  if (prod(g[j:(j + 1)] - target) > 0) next
  lo_hi = ts[j:(j + 1)]
  t = uniroot(function(t) { v = gap(t); if (is.na(v)) NaN else v - target }, lo_hi, tol = 1e-12)$root
  Uy = drop(crossprod(e$vectors, y_of(t))); ex = max(peaks(UX0, Uy, e$values))
  G = cbind(rbinom(n, 2, 0.3)); l = sapply(impl, function(f) f(e$vectors, e$values, Uy, UX0, G, matrix(0, n, 0), FALSE)$l0)
  rows[[length(rows) + 1]] = c(gap = target, l - ex)
}
r = do.call(rbind, rows)
cat(sprintf("%d bimodal cases (%d draws), peak gaps |1e-5 .. 1e-2|\n", nrow(r), tries))
for (k in names(impl)) cat(sprintf("%-17s below exact by > 1e-6: %d   max shortfall %.3g   max |l - baseline| %.3g\n", k,
  sum(r[, k] < -1e-6), -min(r[, k]), max(abs(r[, k] - r[, "baseline"]))))
