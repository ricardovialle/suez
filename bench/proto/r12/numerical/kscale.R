# Kernel-scale invariance (tests/model_checks.R lines 149-154) at kernel level: base12 and prototypes vs an exact
# optimiser over delta = s2k/s2 (x = log10 delta in [-25, 25] + both boundaries) for f = 1, 1e6, 1e-6.
#   Rscript bench/proto/r12/numerical/kscale.R proto.cpp ...
Sys.setenv(SUEZ_HARNESS_LIB = "bench/lib/base12"); source("bench/proto/harness.R")
check = function(...) invisible(NULL); src = readLines("tests/model_checks.R"); eval(parse(text = src[grep("^profile_ml = function", src):grep("^anno = data.frame", src)]))
impl = list(base12 = baseline_fit)
for (f in commandArgs(TRUE)) { e = new.env(); Rcpp::sourceCpp(f, env = e); impl[[sub("\\.cpp$", "", basename(f))]] = e$fit_gene_proto }
exact = function(UX, Uy, lam) {   # sup over delta in [0, inf] (h = 1 only when every lambda > 0)
  n = length(Uy)
  fd = function(d) { W = 1 / d; b = if (ncol(UX)) qr.coef(qr(UX * sqrt(W)), Uy * sqrt(W)) else numeric(0)
    r = Uy - if (ncol(UX)) UX %*% b else 0; -0.5 * n - 0.5 * (n * log(sum(r^2 * W) / n) + sum(log(d))) }
  fx = function(x) fd(10^x * lam + 1)
  xs = seq(-25, 25, length.out = 5001); v = sapply(xs, fx); pk = which(diff(sign(diff(v))) == -2) + 1
  c(v, fd(rep(1, n)), if (all(lam > 0)) fd(lam), vapply(pk, function(i) optimize(fx, xs[c(i - 1, i + 1)], maximum = TRUE, tol = 1e-12)$objective, 0)) |> max()
}
X0 = cbind(std(cond)); C = cbind(cond)
res = list(); ex = list()
for (f in c(1, 1e6, 1e-6)) {
  eg = suez:::kernel_contrast_eigen(f * K); U = eg$vectors; lam = eg$values; UX0 = crossprod(U, X0)
  genes = lapply(1:n_gene, function(g) { Gg = t(G[(g - 1) * n_snp + 1:n_snp, ind]); storage.mode(Gg) = "double"
    list(gene = g, seed = g, Uy = drop(crossprod(U, Y[g, ])), G = Gg[, apply(Gg, 2, sd) > 0, drop = FALSE]) })
  P = list(U = U, lambda = lam, UX0 = UX0, C = C, boot = FALSE, genes = genes)
  for (k in names(impl)) res[[k]][[as.character(f)]] = run_impl(impl[[k]], P)
  ex[[as.character(f)]] = do.call(rbind, lapply(genes, function(g) { Uy = g$Uy
    t(sapply(seq_len(ncol(g$G)), function(s) { x = crossprod(U, std(g$G[, s])); xi = crossprod(U, std(g$G[, s] * cond))
      c(l0 = exact(UX0, Uy, lam), l_geno = exact(cbind(UX0, x), Uy, lam), l_interact = exact(cbind(UX0, x, xi), Uy, lam)) })) }))
}
tab = function(r) do.call(rbind, lapply(r$out, function(o) cbind(l0 = o$l0, o$res[, 1:2, drop = FALSE])))
for (k in names(impl)) {
  L = lapply(res[[k]], tab)
  inv = max(sapply(L[-1], function(x) max(abs(x - L[[1]]))))
  dex = sapply(names(L), function(f) { d = L[[f]] - ex[[f]]; c(min = min(d), max = max(d)) })
  cat(sprintf("%-12s invariance max|dl| %.2g | l - exact (min/max) f=1: %.2g/%.2g  f=1e6: %.2g/%.2g  f=1e-6: %.2g/%.2g\n", k, inv,
              dex[1, 1], dex[2, 1], dex[1, 2], dex[2, 2], dex[1, 3], dex[2, 3]))
}
saveRDS(list(res = res, ex = ex), "bench/proto/r12/numerical/kscale_out.rds")
