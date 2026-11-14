# Kernel-scale and boundary-optimum suite for an installed candidate (added after p_prof2: Newton started from -Inf
# could not return to a grid optimum at h = 0 / h = 1, invisible at f = 1 but 3.7e-6..9.7e-6 below the ML at f = 1e6).
#   Rscript bench/proto/r12/review/kscale12.R <cand lib dir> <ref lib dir> <out.tsv>
# 1 kernel_scale: tests/model_checks.R's data (GRM of 120 donors with 1-2 samples, condition rnorm within donor, genes
#   at h in {0, 0.3, 0.8, 1}), sample_kernel = f K, f = 1, 1e6, 1e-6: l0 / l_geno / l_interact equal across scales
#   within 1e-6 (the test's tolerance) and within 1e-9 of an exact optimiser over log10 delta in [-25, 25] plus both
#   boundaries (delta = 0, and delta = Inf when no eigenvalue is 0).
# 2 boundary_h0 / boundary_h1: genes whose ML is exactly h = 0 (s2k = 0) or h = 1 (s2 = 0; full-rank donor GRM, one sample
#   per donor), kernel scales 1, 1e3, 1e6; l within 1e-9 of the exact optimum.
a = commandArgs(TRUE); Sys.setenv(SUEZ_HARNESS_LIB = a[2]); out = a[3]
source("bench/proto/harness.R")
check = function(...) invisible(NULL); src = readLines("tests/model_checks.R")
eval(parse(text = src[grep("^profile_ml = function", src):grep("^anno = data.frame", src)]))   # nd, ind, K, G, Y, cond, std()
exact = function(UX, Uy, lam) {   # sup over delta = s2k / s2 in [0, Inf] (numerical/kscale.R)
  n = length(Uy)
  fd = function(d) { W = 1 / d; b = if (ncol(UX)) qr.coef(qr(UX * sqrt(W)), Uy * sqrt(W)) else numeric(0)
    r = Uy - if (ncol(UX)) UX %*% b else 0; -0.5 * n - 0.5 * (n * log(sum(r^2 * W) / n) + sum(log(d))) }
  fx = function(x) fd(10^x * lam + 1)
  xs = seq(-25, 25, length.out = 5001); v = sapply(xs, fx); pk = which(diff(sign(diff(v))) == -2) + 1
  b0 = fd(rep(1, n)); b1 = if (all(lam > 0)) fd(lam) else -Inf
  pv = vapply(pk, function(i) optimize(fx, xs[c(i - 1, i + 1)], maximum = TRUE, tol = 1e-12)$objective, 0)
  m = max(v, pv, b0, b1)   # boundary optimum: the sup is attained at delta = 0 or delta = Inf (to 1e-9; the grid ends are plateaus)
  at = if (b0 >= m - 1e-9) 0 else if (b1 >= m - 1e-9) 1 else 0.5
  c(max = m, at = at)
}
X0 = cbind(std(cond)); C = cbind(cond)
# a prepare()-like case: per gene Uy and samples x SNPs genotypes; exact optima of l0, l_geno, l_interact per SNP
mkcase = function(Kf, Ys, Gs, cnd) {
  eg = suez:::kernel_contrast_eigen(Kf); U = eg$vectors; lam = eg$values; UX0 = crossprod(U, cbind(std(cnd)))
  genes = lapply(seq_along(Ys), function(g) list(gene = g, seed = g, Uy = drop(crossprod(U, Ys[[g]])), G = Gs[[g]]))
  ex = do.call(rbind, lapply(genes, function(g) t(sapply(seq_len(ncol(g$G)), function(s) {
    x = crossprod(U, std(g$G[, s])); xi = crossprod(U, std(g$G[, s] * cnd))
    c(l0 = exact(UX0, g$Uy, lam)[["max"]], l_geno = exact(cbind(UX0, x), g$Uy, lam)[["max"]],
      l_interact = exact(cbind(UX0, x, xi), g$Uy, lam)[["max"]]) }))))
  list(P = list(U = U, lambda = lam, UX0 = UX0, C = cbind(cnd), boot = FALSE, genes = genes), ex = ex, U = U, lam = lam, UX0 = UX0)
}
cases = list()
Gm = lapply(1:n_gene, function(g) { x = t(G[(g - 1) * n_snp + 1:n_snp, ind]); storage.mode(x) = "double"; x[, apply(x, 2, sd) > 0, drop = FALSE] })
for (f in c(1, 1e6, 1e-6)) cases[[sprintf("kernel_scale f=%g", f)]] = mkcase(f * K, lapply(1:n_gene, function(g) Y[g, ]), Gm, cond)
# boundary optima: draw y in the contrast basis until the exact l0 optimum is the wanted boundary
set.seed(11)
K1 = grm[inds, inds]; cond1 = rnorm(nd)                          # one sample per donor: no zero contrast eigenvalue
G1 = lapply(1:4, function(g) { x = matrix(rbinom(nd * 3, 2, 0.3), nd); storage.mode(x) = "double"; x[, apply(x, 2, sd) > 0, drop = FALSE] })
draw = function(Kf, cnd, at, n_genes = 4) {
  e = suez:::kernel_contrast_eigen(Kf); UX0 = crossprod(e$vectors, cbind(std(cnd))); ys = list(); tries = 0
  while (length(ys) < n_genes && tries < 500) {
    tries = tries + 1; lam = e$values / mean(e$values)
    z = rnorm(length(lam)) * if (at == 0) ifelse(lam > median(lam), 0.6, 1) else sqrt(lam) * ifelse(lam < median(lam), 0.6, 1)
    y = drop(e$vectors %*% z) + 0.3 * std(cnd)
    if (exact(UX0, drop(crossprod(e$vectors, y)), e$values)[["at"]] == at) ys[[length(ys) + 1]] = y
  }
  ys
}
for (at in c(0, 1)) {
  Kb = if (at == 1) K1 else K; cb = if (at == 1) cond1 else cond
  ys = draw(Kb, cb, at); if (!length(ys)) stop("no boundary case h = ", at)
  Gb = if (at == 1) G1 else Gm[1:4]
  for (f in c(1, 1e3, 1e6)) cases[[sprintf("boundary_h%d f=%g", at, f)]] = mkcase(f * Kb, ys, Gb[seq_along(ys)], cb)
}
fi = tempfile(fileext = ".rds"); fo = tempfile(fileext = ".rds"); saveRDS(lapply(cases, `[[`, "P"), fi)
system2("Rscript", c("bench/proto/review/run_lib.R", a[1], fi, fo)); cr = readRDS(fo)
tab = function(r) do.call(rbind, lapply(r$out, function(o) cbind(l0 = o$l0, o$res[, 1:2, drop = FALSE])))
rows = list(); L = list()
for (nm in names(cases)) {
  rb = run_impl(baseline_fit, cases[[nm]]$P); rc = cr[[nm]]
  if (inherits(rc, "error")) { rows[[nm]] = data.frame(dataset = nm, check = "run", value = NA, threshold = "-", status = "FAIL", detail = conditionMessage(rc)); next }
  A = tab(rb); B = tab(rc); E = cases[[nm]]$ex; L[[nm]] = B
  db = max(abs(A - E)); dc = max(abs(B - E)); na = sum(is.na(A) != is.na(B))
  gate = !startsWith(nm, "kernel_scale") || grepl("f=1e\\+06", nm)   # the test's f = 1 / 1e-6 exact gaps: reported only
  st = if (na || is.na(dc)) "FAIL" else if (dc <= 1e-9) "PASS" else if (gate) "FAIL" else if (dc <= db + 1e-9) "PASS" else if (dc <= 1e-6) "WARN" else "FAIL"   # f = 1e-6: base12 itself is 9e-8 off
  rows[[nm]] = data.frame(dataset = nm, check = "max |l - exact ML| (l0, l_geno, l_interact)", value = signif(dc, 4),
    threshold = if (gate) "1e-9" else "1e-9 or <= ref + 1e-9 (WARN to 1e-6)", status = st,
    detail = sprintf("ref %.3g; cand min/max l - exact %.3g/%.3g; %d fits; NA mismatches %d", db, min(B - E), max(B - E), length(E), na))
}
inv = max(sapply(c("kernel_scale f=1e+06", "kernel_scale f=1e-06"), function(k) max(abs(L[[k]] - L[["kernel_scale f=1"]]))))
rows[["inv"]] = data.frame(dataset = "kernel_scale", check = "max |dl| across f = 1, 1e6, 1e-6 (tests/model_checks.R)", value = signif(inv, 4),
                           threshold = "1e-6", status = if (is.finite(inv) && inv < 1e-6) "PASS" else "FAIL", detail = "")
write.table(do.call(rbind, rows), out, sep = "\t", quote = FALSE, row.names = FALSE)
