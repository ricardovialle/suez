# Edge cases of tests/model_checks.R run through a prototype and the baseline kernel:
#  A. centred donor GRM (K 1 = 0) with repeated samples, condition varying within individuals, genes with h in
#     {0, 0.3, 0.8, 1} (boundaries s2k = 0 and s2 ~ 0), bootstrap on;
#  B. bimodal likelihood in h (two near-equal peaks, 2e-3 apart), checked against a dense exact profile;
#  C. same as A but y with no within-individual variation at all (ML at the s2 -> 0 end, Newton up to x = 12).
# Rscript bench/proto/algorithm/edge_checks.R <proto.cpp>
a = commandArgs(TRUE); source("bench/proto/harness.R")
Rcpp::sourceCpp(a[1], cacheDir = "bench/proto/algorithm/.cache")
std = function(x) (x - mean(x)) / sd(x)
cmp = function(name, A, B) {
  d = sapply(c("l0"), function(k) max(abs(sapply(A, `[[`, k) - sapply(B, `[[`, k))))
  RA = do.call(rbind, lapply(A, `[[`, "res")); RB = do.call(rbind, lapply(B, `[[`, "res"))
  dd = apply(abs(RA - RB), 2, max, na.rm = TRUE)
  st = if (exists("proto_stats")) paste(proto_stats(TRUE), collapse = "/") else ""
  cat(sprintf("[nonNA l_geno %d/%d, l0 %.4g, stats %s] ", sum(!is.na(RA[, 1])), nrow(RA), A[[1]]$l0, st))
  cat(sprintf("%-28s l0 %.2g | l_geno %.2g l_int %.2g beta %.2g se %.2g lb_geno %.2g lb_int %.2g | NA mismatch %d\n", name, d,
              dd[1], dd[2], dd[4], dd[5], dd[6], dd[7], sum(is.na(RA) != is.na(RB))))
}
runboth = function(U, lam, Ys, X0, Gl, C, boot) {
  UX0 = crossprod(U, X0); seeds = 1:length(Ys)
  f = function(fun) lapply(seq_along(Ys), function(i) { set.seed(seeds[i]); fun(U, lam, drop(crossprod(U, Ys[[i]])), UX0, Gl[[i]], C, boot) })
  list(base = f(baseline_fit), proto = f(fit_gene_proto))
}
# --- A / C
set.seed(7)
nd = 120; reps = sample(1:2, nd, TRUE, prob = c(0.5, 0.5))
inds = sprintf("I%03d", 1:nd); ind = rep(inds, reps); N = length(ind)
Am = matrix(rbinom(nd * 2000, 2, 0.3), nd)
Z = scale(Am); grm = tcrossprod(Z) / ncol(Z); dimnames(grm) = list(inds, inds)
K = grm[ind, ind]
cond = rnorm(N); n_gene = 8; n_snp = 6
G = matrix(rbinom(n_gene * n_snp * nd, 2, 0.3), n_gene * n_snp, nd, dimnames = list(paste0("rs", 1:(n_gene * n_snp)), inds))
u = drop(t(chol(K + 1e-6 * diag(N))) %*% rnorm(N))
Y = lapply(1:n_gene, function(g) {
  h = c(0, 0.3, 0.8, 1)[g %% 4 + 1]; gs = G[(g - 1) * n_snp + 1, ind]
  0.3 * cond + 0.4 * gs + 0.3 * gs * cond + sqrt(h) * std(u) + sqrt(1 - h) * rnorm(N) })
eig = suez:::kernel_contrast_eigen(K)
Gl = lapply(1:n_gene, function(g) { x = t(G[(g - 1) * n_snp + 1:n_snp, ind]); storage.mode(x) = "double"; x })
r = runboth(eig$vectors, eig$values, Y, cbind(std(cond)), Gl, cbind(cond), TRUE); cmp("A: GRM, h in {0,.3,.8,1}", r$base, r$proto)
# C: y constant within individual (random effect + fixed effects of individual-level terms only)
uind = rnorm(nd); names(uind) = inds; gi = lapply(1:n_gene, function(g) G[(g - 1) * n_snp + 1, ind])
Yc = lapply(1:n_gene, function(g) 0.4 * gi[[g]] + uind[ind] + 1e-9 * (g %% 2) * rnorm(N))
condi = rnorm(nd); names(condi) = inds; ci = condi[ind]
r = runboth(eig$vectors, eig$values, Yc, cbind(std(ci)), Gl, cbind(ci), TRUE); cmp("C: no within-indiv. noise", r$base, r$proto)
Kd = outer(ind, ind, "=="); class(Kd) = "numeric"; eigd = suez:::kernel_contrast_eigen(Kd)
r = runboth(eigd$vectors, eigd$values, Yc, cbind(std(ci)), Gl, cbind(ci), TRUE); cmp("C': same, default kernel", r$base, r$proto)
r = runboth(eigd$vectors, eigd$values, Y, cbind(std(cond)), Gl, cbind(cond), TRUE); cmp("A': A with default kernel", r$base, r$proto)
# --- B: bimodal
exact_ml = function(UX, Uy, lam) {
  f = function(x) { h = 1 / (1 + 10^-x); d = h * lam + (1 - h); W = 1 / d
    b = solve(crossprod(UX, UX * W), crossprod(UX, Uy * W))
    -0.5 * length(Uy) - 0.5 * sum(log(sum((Uy - UX %*% b)^2 * W) / length(Uy) * d)) }
  xs = seq(-8, 8, length.out = 1601); v = sapply(xs, f); pk = which(diff(sign(diff(v))) == -2) + 1
  sapply(pk, function(i) optimize(f, xs[c(i - 1, i + 1)], maximum = TRUE, tol = 1e-10)$objective)
}
set.seed(1859)
nb = 120; Qb = qr.Q(qr(matrix(1, nb, 1)), complete = TRUE)[, -1]
lb = c(rep(1e-3, 40), rep(1, 40), rep(30, 39)) * exp(rnorm(nb - 1, 0, 0.3))
vb = c(runif(1, 0.01, 2), runif(1, 0.01, 2), runif(1, 0.01, 60))
zb = rnorm(nb - 1) * sqrt(rep(vb, c(40, 40, 39)))
Kb = Qb %*% (lb * t(Qb)) + 1 / nb
eb = suez:::kernel_contrast_eigen(Kb); condb = rnorm(nb); UX0b = crossprod(eb$vectors, std(condb))
y_of = function(t) drop(Qb %*% (zb * rep(c(1, 1, t), c(40, 40, 39))))
peaks = function(t) exact_ml(UX0b, drop(crossprod(eb$vectors, y_of(t))), eb$values)
sgap = function(t) { p = peaks(t); p[1] - p[length(p)] }
tb = sapply(c(2e-3, -2e-3), function(g) uniroot(function(t) sgap(t) - g, c(0.9, 1))$root)
Gb = matrix(rbinom(4 * nb, 2, 0.3), nb, 4); storage.mode(Gb) = "double"
Yb = list(y_of(tb[1]), y_of(tb[2]))
r = runboth(eb$vectors, eb$values, Yb, cbind(std(condb)), list(Gb, Gb), matrix(0, nb, 0), FALSE); cmp("B: bimodal (eQTL, no GxE)", r$base, r$proto)
db = sapply(1:2, function(k) { Uy = drop(crossprod(eb$vectors, Yb[[k]]))
  p0 = exact_ml(UX0b, Uy, eb$values)
  c(length(p0), r$proto[[k]]$l0 - max(p0), sapply(1:4, function(j) r$proto[[k]]$res[j, 1] - max(exact_ml(cbind(UX0b, crossprod(eb$vectors, std(Gb[, j]))), Uy, eb$values)))) })
cat("B vs exact profile: peaks", db[1, ], " max |dl|", max(abs(db[-1, ])), "\n")
