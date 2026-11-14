# Global-optimum stress for an installed candidate (stress_exact.R / stress_bimodal.R recipes, out of process):
# random 3-block kernel spectra (n = 120; often bimodal in h, narrow peaks, boundary optima) and spectra tuned so
# two peaks differ by 1e-5..1e-2. l0 and l_geno of SNP 1 vs the exact ML (1601-point grid in log10 delta, every
# local maximum refined); reference and candidate must not miss more often, and agree within 1e-9 (WARN 1e-6).
#   Rscript bench/proto/r12/review/stress12.R <cand lib dir> <ref lib dir> <out.tsv> [n_random=150] [n_bimodal=20]
a = commandArgs(TRUE); Sys.setenv(SUEZ_HARNESS_LIB = a[2]); out = a[3]
nr = if (length(a) >= 4) as.integer(a[4]) else 150; nb = if (length(a) >= 5) as.integer(a[5]) else 20
source("bench/proto/harness.R")
f_ll = function(UX, Uy, lam) function(x) { h = 1 / (1 + 10^-x); d = h * lam + (1 - h); W = 1 / d
  b = solve(crossprod(UX, UX * W), crossprod(UX, Uy * W)); -0.5 * length(Uy) - 0.5 * sum(log(sum((Uy - UX %*% b)^2 * W) / length(Uy) * d)) }
exact = function(UX, Uy, lam, npts = 1601) {
  f = f_ll(UX, Uy, lam); xs = seq(-8, 8, length.out = npts); v = sapply(xs, f); pk = which(diff(sign(diff(v))) == -2) + 1
  max(c(v[c(1, npts)], vapply(pk, function(i) optimize(f, xs[c(i - 1, i + 1)], maximum = TRUE, tol = 1e-11)$objective, 0)))
}
n = 120; Q = qr.Q(qr(matrix(1, n, 1)), complete = TRUE)[, -1]
mkcase = function(lam, y, kind, cond = rnorm(n)) {
  K = Q %*% (lam * t(Q)) + 1 / n; e = suez:::kernel_contrast_eigen(K)
  UX0 = crossprod(e$vectors, (cond - mean(cond)) / sd(cond)); Uy = drop(crossprod(e$vectors, y))
  G = cbind(rbinom(n, 2, 0.3), rbinom(n, 2, 0.4)); g1 = crossprod(e$vectors, (G[, 1] - mean(G[, 1])) / sd(G[, 1]))
  list(U = e$vectors, lambda = e$values, UX0 = UX0, C = matrix(0, n, 0), boot = FALSE, kind = kind,
       genes = list(list(gene = "g", seed = 1, Uy = Uy, G = G)),
       ex = c(exact(UX0, Uy, e$values), exact(cbind(UX0, g1), Uy, e$values)))
}
set.seed(99); cases = list()
for (i in seq_len(nr)) {   # stress_exact.R
  m = sort(sample(10:90, 2)); sz = c(m[1], m[2] - m[1], n - 1 - m[2])
  lam = rep(10^runif(3, -4, 3), sz) * exp(rnorm(n - 1, 0, 0.3)); z = rnorm(n - 1) * sqrt(rep(10^runif(3, -2, 2), sz))
  cases[[length(cases) + 1]] = mkcase(lam, drop(Q %*% z), "random 3-block spectrum")
}
set.seed(17); sz = c(40, 40, 39); tries = 0
while (sum(sapply(cases, `[[`, "kind") == "bimodal, tuned peak gap") < nb && tries < 30 * nb) {   # stress_bimodal.R
  tries = tries + 1
  lam = rep(10^c(runif(1, -3.5, -2.5), runif(1, -0.5, 0.5), runif(1, 1, 2)), sz) * exp(rnorm(n - 1, 0, 0.3))
  K = Q %*% (lam * t(Q)) + 1 / n; e = suez:::kernel_contrast_eigen(K); cond = rnorm(n)
  UX0 = crossprod(e$vectors, (cond - mean(cond)) / sd(cond))
  z = rnorm(n - 1) * sqrt(rep(c(runif(2, 0.01, 2), runif(1, 0.01, 60)), sz)); y_of = function(t) drop(Q %*% (z * rep(c(1, 1, t), sz)))
  peaks = function(Uy) { f = f_ll(UX0, Uy, e$values); xs = seq(-8, 8, length.out = 801); v = sapply(xs, f); pk = which(diff(sign(diff(v))) == -2) + 1
    vapply(pk, function(i) optimize(f, xs[c(i - 1, i + 1)], maximum = TRUE, tol = 1e-11)$objective, 0) }
  gap = function(t) { p = peaks(drop(crossprod(e$vectors, y_of(t)))); if (length(p) < 2) NA else p[1] - p[length(p)] }
  ts = 10^seq(-0.5, 0.5, length.out = 9); g = sapply(ts, gap)
  j = which(!is.na(g[-1]) & !is.na(g[-9]) & sign(g[-1]) != sign(g[-9]))[1]; if (is.na(j)) next
  target = sample(c(-1, 1), 1) * 10^runif(1, -5, -2); if (prod(g[j:(j + 1)] - target) > 0) next
  t = tryCatch(uniroot(function(t) { v = gap(t); if (is.na(v)) NaN else v - target }, ts[j:(j + 1)], tol = 1e-12)$root, error = function(e) NA)
  if (is.na(t)) next
  cases[[length(cases) + 1]] = mkcase(lam, y_of(t), "bimodal, tuned peak gap", cond)   # same cond: the tuned peaks are those of l0
}
fi = tempfile(fileext = ".rds"); fo = tempfile(fileext = ".rds"); saveRDS(cases, fi)
system2("Rscript", c("bench/proto/review/run_lib.R", a[1], fi, fo)); cr = readRDS(fo)
res = do.call(rbind, lapply(seq_along(cases), function(i) {
  P = cases[[i]]; b = run_impl(baseline_fit, P)$out[[1]]; c = cr[[i]]$out[[1]]
  data.frame(kind = P$kind, base_l0 = b$l0 - P$ex[1], cand_l0 = c$l0 - P$ex[1], base_l1 = b$res[1, 1] - P$ex[2], cand_l1 = c$res[1, 1] - P$ex[2],
             d = max(abs(c(c$l0 - b$l0, c$res[, 1] - b$res[, 1])), na.rm = TRUE), na = sum(is.na(c$res[, 1]) != is.na(b$res[, 1])))
}))
rows = lapply(split(res, res$kind), function(r) {
  mb = sum(r$base_l0 < -1e-6, r$base_l1 < -1e-6, na.rm = TRUE); mc = sum(r$cand_l0 < -1e-6, r$cand_l1 < -1e-6, na.rm = TRUE)
  st = if (mc > mb || max(r$d) > 1e-6 || sum(r$na)) "FAIL" else if (max(r$d) > 1e-9) "WARN" else "PASS"
  data.frame(dataset = r$kind[1], check = "fits below the exact ML by > 1e-6 (cand)", value = mc, threshold = "<= ref; |cand-ref| <= 1e-9",
             status = st, detail = sprintf("%d cases; ref misses %d; max shortfall ref %.3g cand %.3g; max |cand-ref| %.3g; NA mismatches %d",
                                           nrow(r), mb, -min(c(r$base_l0, r$base_l1), na.rm = TRUE), -min(c(r$cand_l0, r$cand_l1), na.rm = TRUE), max(r$d), sum(r$na)))
})
write.table(do.call(rbind, rows), out, sep = "\t", quote = FALSE, row.names = FALSE)
