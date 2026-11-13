# Near-collinear but independent columns: cov2 = cov1 + eps * noise (X0 block) and a SNP column = covariate +
# eps * noise. 1 - R^2 ~ eps^2, so p9's COLLINEAR = 1e-12 pivot rule acts at eps <~ 1e-6. l0 / l_geno of the
# baseline and a candidate vs exact ML with a full-rank (LAPACK, no rank cut) QR.
#   Rscript bench/proto/review/cand/collinear_sweep.R cand.cpp
source("bench/proto/harness.R"); e = new.env(); Rcpp::sourceCpp(commandArgs(TRUE)[1], env = e); cand = e$fit_gene_proto
exact = function(UX, Uy, lam) {
  f = function(x) { h = 1 / (1 + 10^-x); d = h * lam + (1 - h); W = 1 / d
    b = qr.coef(qr(UX * sqrt(W), LAPACK = TRUE), Uy * sqrt(W)); r = Uy - UX %*% b
    -0.5 * length(Uy) - 0.5 * sum(log(sum(r^2 * W) / length(Uy) * d)) }
  xs = seq(-8, 8, length.out = 801); v = sapply(xs, f); i = which.max(v)
  max(v[i], optimize(f, xs[c(max(1, i - 1), min(801, i + 1))], maximum = TRUE, tol = 1e-12)$objective)
}
set.seed(8); ni = 150; ind = rep(1:ni, each = 2); N = 2 * ni; cnd = rep(c(0, 1), ni)
Z = matrix(rbinom(ni * 2000, 2, 0.3), ni); Z = scale(Z[, apply(Z, 2, sd) > 0]); K = (tcrossprod(Z) / ncol(Z))[ind, ind]
eig = suez:::kernel_contrast_eigen(K); U = eig$vectors; lam = eig$values
g = rbinom(ni, 2, 0.3)[ind]; c1 = rnorm(N); y = 0.3 * cnd + 0.3 * c1 + 0.3 * g + rnorm(ni)[ind] + rnorm(N)
std = function(x) (x - mean(x)) / sd(x)
cat(sprintf("%-6s %-28s %10s %10s %10s\n", "eps", "case", "base-exact", "cand-exact", "cand-base"))
for (eps in 10^-(3:9)) {
  X0a = cbind(std(cnd), std(c1), std(c1 + eps * rnorm(N)))                   # near-duplicate covariates
  X0b = cbind(std(cnd), std(g + eps * rnorm(N)))                              # covariate ~ the SNP
  for (cs in list(list("cov2 = cov1 + eps", X0a), list("covariate = SNP + eps", X0b))) {
    UX0 = crossprod(U, cs[[2]]); Uy = drop(crossprod(U, y)); G = cbind(g)
    b = baseline_fit(U, lam, Uy, UX0, G, cbind(cnd), FALSE); p = cand(U, lam, Uy, UX0, G, cbind(cnd), FALSE)
    e0 = exact(UX0, Uy, lam); e1 = exact(cbind(UX0, crossprod(U, std(g))), Uy, lam)
    cat(sprintf("%-6.0e %-28s %10.2g %10.2g %10.2g   (l_geno: %.2g %.2g)\n", eps, cs[[1]], b$l0 - e0, p$l0 - e0, p$l0 - b$l0,
                b$res[1, 1] - e1, p$res[1, 1] - e1))
  }
}
