# Extreme cases for eigenvalue-grouping candidates: null directions with h -> 1, kernels scaled x1e6 / x1e-6,
# null eigenvalues spread by rounding noise. Each impl vs the baseline and vs an exact optimizer on the true
# per-contrast lambda (x in [-12, 12], 4801 points, every local max refined, plus h = 0 and h = 1).
#   Rscript bench/proto/review/cand/extreme_cases.R cand1.cpp [cand2.cpp ...]
source("bench/proto/harness.R")
files = commandArgs(TRUE)
impl = list(baseline = baseline_fit)
for (f in files) { e = new.env(); Rcpp::sourceCpp(f, env = e); impl[[sub("\\.cpp$", "", basename(f))]] = e$fit_gene_proto }
std = function(x) (x - mean(x)) / sd(x)
mk = function(K, cond, G, ys) {
  eig = suez:::kernel_contrast_eigen(K); X0 = suez:::standardize_cols(model.matrix(~ c, data.frame(c = cond))[, -1, drop = FALSE])
  C = model.matrix(~ c, data.frame(c = cond))[, -1, drop = FALSE]
  list(U = eig$vectors, lambda = eig$values, UX0 = crossprod(eig$vectors, X0), C = C, boot = TRUE,
       genes = lapply(seq_along(ys), function(i) list(gene = paste0("g", i), seed = 10 + i, Uy = drop(crossprod(eig$vectors, ys[[i]])), G = G)))
}
exact = function(UX, Uy, lam) {
  f = function(h) { d = h * lam + (1 - h); if (any(d <= 0)) return(-Inf); W = 1 / d
    b = if (ncol(UX)) qr.coef(qr(UX * sqrt(W)), Uy * sqrt(W)) else numeric(0); r = Uy - if (ncol(UX)) UX %*% b else 0
    -0.5 * length(Uy) - 0.5 * sum(log(sum(r^2 * W) / length(Uy) * d)) }
  fx = function(x) f(1 / (1 + 10^-x)); xs = seq(-12, 12, length.out = 4801); v = sapply(xs, fx)
  pk = which(diff(sign(diff(v))) == -2) + 1
  max(c(f(0), f(1), v, vapply(pk, function(i) optimize(fx, xs[c(i - 1, i + 1)], maximum = TRUE, tol = 1e-12)$objective, 0)))
}
set.seed(5)
ni = 100; ind = rep(sprintf("i%03d", 1:ni), each = 2); N = 2 * ni; cnd = rep(c("a", "b"), ni)
Kd = outer(ind, ind, "==") * 1
Z = matrix(rbinom(ni * 2000, 2, 0.3), ni); Z = scale(Z[, apply(Z, 2, sd) > 0]); Kg = tcrossprod(Z) / ncol(Z); Kg = Kg[rep(1:ni, each = 2), rep(1:ni, each = 2)]
G = sapply(runif(8, 0.1, 0.5), function(f) rbinom(ni, 2, f))[rep(1:ni, each = 2), ]
u = rnorm(ni)[rep(1:ni, each = 2)]
ys = function(within) lapply(within, function(e) 0.3 * (cnd == "b") + u + 0.2 * G[, 1] + e * rnorm(N))
noisy = function(K, eps) { E = matrix(rnorm(N * N, 0, eps), N); K + (E + t(E)) / 2 }
cases = list(
  default_x1e6_h_to_1     = mk(Kd * 1e6, cnd, G, ys(c(1, 1e-2, 1e-4))),
  default_x1em6_h_to_1    = mk(Kd * 1e-6, cnd, G, ys(c(1, 1e-2, 1e-4))),
  default_h_to_1          = mk(Kd, cnd, G, ys(c(1e-2, 1e-4, 1e-6))),
  grm_reps_h_to_1         = mk(Kg, cnd, G, ys(c(1e-2, 1e-4, 1e-6))),
  grm_reps_x1e6_h_to_1    = mk(Kg * 1e6, cnd, G, ys(c(1, 1e-3, 1e-5))),
  null_noise_1em13        = mk(noisy(Kg, 1e-13), cnd, G, ys(c(1, 1e-3, 1e-5))),
  null_noise_1em11        = mk(noisy(Kg, 1e-11), cnd, G, ys(c(1, 1e-3, 1e-5))),
  default_x1e6_opt_x4     = mk(Kd * 1e6, cnd, G, ys(c(3e-5, 1e-5, 3e-6, 1e-6))),
  grm_reps_x1e6_opt_x4    = mk(Kg * 1e6, cnd, G, ys(c(3e-5, 1e-5, 3e-6, 1e-6))),
  default_x1e3_opt_x4     = mk(Kd * 1e3, cnd, G, ys(c(1e-3, 3e-4, 1e-4))),
  near_perfect_fit        = mk(Kd, cnd, G, lapply(c(1e-3, 1e-6, 1e-9), function(e) (cnd == "b") + e * rnorm(N))))
if (nzchar(Sys.getenv("ONLY"))) cases = cases[strsplit(Sys.getenv("ONLY"), ",")[[1]]]
options(width = 200)
for (nm in names(cases)) {
  P = cases[[nm]]; cat(sprintf("\n## %s  (null eigenvalues: %d exact 0, %d in (0, 1e-8 max], max %.2g)\n", nm, sum(P$lambda == 0),
                               sum(P$lambda > 0 & P$lambda <= 1e-8 * max(P$lambda)), max(P$lambda[P$lambda <= 1e-8 * max(P$lambda)])))
  R = lapply(impl, function(f) tryCatch(run_impl(f, P), error = function(e) e))
  ex = sapply(P$genes, function(g) c(l0 = exact(P$UX0, g$Uy, P$lambda), l1 = exact(cbind(P$UX0, crossprod(P$U, std(g$G[, 1]))), g$Uy, P$lambda)))
  tab = t(sapply(names(R), function(k) { r = R[[k]]; if (inherits(r, "error")) return(c(rep(NA, 5)))
    l0 = sapply(r$out, `[[`, "l0"); l1 = sapply(r$out, function(o) o$res[1, 1])
    cmp = if (k == "baseline") rep(0, 3) else compare(R$baseline, r)[c("l_geno", "l_boot_geno", "na_mismatch")]
    c(max_l0_minus_exact = max(abs(l0 - ex[1, ])), max_l1_minus_exact = max(abs(l1 - ex[2, ])), vs_base_l_geno = cmp[[1]], vs_base_l_boot = cmp[[2]], na_mismatch = cmp[[3]]) }))
  print(signif(tab, 3))
}
