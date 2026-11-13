# Edge-case suite for candidate fitters: runs a prototype fit_gene_proto (Rcpp::sourceCpp file with the
# signature of bench/proto/p0_baseline.cpp) and the installed baseline .Call on inputs built to hit the
# numerically fragile parts of src/lmm.cpp, and reports compare() per case.
#   Rscript bench/proto/review/edge_cases.R [candidate.cpp]     (default: bench/proto/p0_baseline.cpp)
# PASS = CONTEXT.md tolerances (|dl| <= 1e-6, beta/se rel <= 1e-6, same df and NA pattern).
# Cases tagged "report" have no well-defined answer (rank-deficient designs, near-perfect fits): a
# difference there is not a bug by itself, but must be listed in the candidate's report.
source("bench/proto/harness.R")
a = commandArgs(TRUE); cand = if (length(a)) a[1] else "bench/proto/p0_baseline.cpp"
# "lib:<dir>" = an installed build, run out of process by run_lib.R (its results are fetched per case below)
if (startsWith(cand, "lib:")) fit_gene_proto = NULL else Rcpp::sourceCpp(cand)

std = function(x) (x - mean(x)) / sd(x)
# a prepare()-like list from a kernel, sample annotation, samples x SNPs genotypes and expression vectors
mk = function(K, ind, cond, G, ys, covs = NULL, boot = TRUE, raw_X0 = FALSE) {
  eig = suez:::kernel_contrast_eigen(K)
  one = length(unique(cond)) == 1
  C = if (one) matrix(0, length(ind), 0) else if (is.numeric(cond)) cbind(cond) else
    model.matrix(~ condition, data.frame(condition = cond))[, -1, drop = FALSE]
  X0 = if (one) matrix(0, length(ind), 0) else model.matrix(~ condition, data.frame(condition = cond))
  if (!is.null(covs)) X0 = cbind(X0, covs)
  X0 = X0[, apply(X0, 2, sd) > 0, drop = FALSE]
  if (!raw_X0) X0 = suez:::standardize_cols(X0)
  if (!is.list(G)) G = rep(list(G), length(ys))
  genes = lapply(seq_along(ys), function(i) list(gene = paste0("g", i), seed = 1000 + i,
                 Uy = drop(crossprod(eig$vectors, ys[[i]])), G = G[[i]]))
  list(U = eig$vectors, lambda = eig$values, UX0 = crossprod(eig$vectors, X0), C = C, boot = boot, genes = genes)
}
geno = function(n, S, af = runif(S, 0.05, 0.5)) sapply(af, function(f) rbinom(n, 2, f))
ind_kernel = function(ind) { K = outer(ind, ind, "=="); K[] = as.numeric(K); K }
grm = function(ni, ind) { Z = geno(ni, 2000); Z = scale(Z[, apply(Z, 2, sd) > 0]); K = tcrossprod(Z) / ncol(Z); dimnames(K) = list(unique(ind), unique(ind)); K[ind, ind] }

set.seed(2024)
cases = list(); tag = list()
add = function(name, P, kind = "exact") { cases[[name]] <<- P; tag[[name]] <<- kind }

# 1 default kernel, balanced 2 reps, condition within individual: eigenvalues exactly {0, 2}, h -> 1 hits d = 0
ni = 100; ind = rep(sprintf("i%03d", 1:ni), each = 2); N = length(ind); cnd = rep(c("a", "b"), ni)
Gi = geno(ni, 40); Gs = Gi[match(ind, unique(ind)), ]
u = rnorm(ni)[match(ind, unique(ind))]
ys = lapply(1:6, function(k) 0.3 * (cnd == "b") + c(0, 0.3, 1, 3, 10, 0)[k] * u + 0.4 * Gs[, 1] * (k == 6) + rnorm(N))
add("default_kernel_balanced", mk(ind_kernel(ind), ind, cnd, Gs, ys))
# 2 unbalanced reps (1-4 per individual): many distinct eigenvalues, few null directions
reps = sample(1:4, ni, TRUE); indu = rep(sprintf("i%03d", 1:ni), reps); Nu = length(indu)
cndu = sample(c("a", "b"), Nu, TRUE); Gu = Gi[match(indu, unique(ind)), ]
uu = rnorm(ni)[match(indu, unique(ind))]
add("default_kernel_unbalanced", mk(ind_kernel(indu), indu, cndu, Gu, lapply(1:4, function(k) k * 0.3 * uu + rnorm(Nu))))
# 3 centred GRM (K 1 = 0), numeric condition varying within individuals
Kg = grm(ni, ind); cn = rnorm(N); ug = drop(t(chol(Kg + 1e-6 * diag(N))) %*% rnorm(N))
ysg = lapply(c(0, 0.5, 0.9, 0.999), function(h) sqrt(h) * std(ug) + sqrt(1 - h) * rnorm(N) + 0.3 * cn + 0.3 * Gs[, 2] * cn)
add("grm_numeric_cond_h0_to_1", mk(Kg, ind, cn, Gs, ysg))
# 4 kernel eigenvalues over 8 orders of magnitude, and the same kernel x 1e6 and x 1e-6 (grid edges |x| = 5)
Q = qr.Q(qr(matrix(rnorm(N * N), N))); ev = 10^seq(-4, 4, length.out = N); Kw = Q %*% (ev * t(Q)); Kw = (Kw + t(Kw)) / 2
yw = lapply(1:4, function(k) drop(Q %*% (sqrt(ev * c(0.01, 0.1, 1, 10)[k] + 1) * rnorm(N))))
for (s in c(1, 1e6, 1e-6)) add(sprintf("wide_spectrum_kernel_x%g", s), mk(Kw * s, ind, cnd, Gs, yw))
# 5 bimodal likelihood in h (construction of tests/model_checks.R: two near-equal peaks)
nb = 120; Qb = qr.Q(qr(matrix(1, nb, 1)), complete = TRUE)[, -1]
lb = c(rep(1e-3, 40), rep(1, 40), rep(30, 39)) * exp(rnorm(nb - 1, 0, 0.3)); Kb = Qb %*% (lb * t(Qb)) + 1 / nb
eb = suez:::kernel_contrast_eigen(Kb); condb = rnorm(nb); UXb = crossprod(eb$vectors, std(condb))
exact_peaks = function(UX, Uy, lam, ends = FALSE) {   # every local maximum in x = log10(delta), refined (tests/model_checks.R)
  f = function(x) { h = 1 / (1 + 10^-x); d = h * lam + (1 - h); W = 1 / d
    b = solve(crossprod(UX, UX * W), crossprod(UX, Uy * W))
    -0.5 * length(Uy) - 0.5 * sum(log(sum((Uy - UX %*% b)^2 * W) / length(Uy) * d)) }
  xs = seq(-8, 8, length.out = 1601); v = sapply(xs, f); pk = which(diff(sign(diff(v))) == -2) + 1
  c(vapply(pk, function(i) optimize(f, xs[c(i - 1, i + 1)], maximum = TRUE, tol = 1e-10)$objective, 0), if (ends) v[c(1, length(v))])
}
repeat {   # variance profile with two peaks, then tune the third block's scale t so the peaks differ by `gap`
  zb = rnorm(nb - 1) * sqrt(rep(c(runif(2, 0.01, 2), runif(1, 0.01, 60)), c(40, 40, 39)))
  y_of = function(t) drop(Qb %*% (zb * rep(c(1, 1, t), c(40, 40, 39))))
  sgap = function(t) { p = exact_peaks(UXb, drop(crossprod(eb$vectors, y_of(t))), eb$values); if (length(p) < 2) NA else p[1] - p[length(p)] }
  g = sapply(c(0.9, 1), sgap); if (!anyNA(g) && prod(g - 2e-3) < 0 && prod(g + 2e-3) < 0) break
}
tb = sapply(c(1e-2, 2e-3, 1e-4, -1e-4, -2e-3, -1e-2), function(g) uniroot(function(t) sgap(t) - g, c(0.9, 1))$root)
add("bimodal_likelihood", mk(Kb, sprintf("b%03d", 1:nb), condb, geno(nb, 30), lapply(tb, y_of)))
# 6 flat likelihood: K = I (default kernel, 1 sample per individual), every grid point ties
ind1 = sprintf("i%03d", 1:150); c1 = sample(c("a", "b"), 150, TRUE)
add("flat_likelihood_identity_kernel", mk(ind_kernel(ind1), ind1, c1, geno(150, 20), lapply(1:3, function(k) rnorm(150))))
# 7 rare SNPs: exactly 5 ALT alleles, carriers all in one condition (interaction column 0 -> df 0), factor3
#   with df 0/1/2 mixed inside one gene
cnd3 = rep(c("a", "b", "c"), length.out = N)
Gr = matrix(0, N, 12)
for (j in 1:12) { pool = which(cnd3 == c("a", "b", "c", "a")[(j - 1) %% 4 + 1]); if (j > 8) pool = seq_len(N)
  Gr[sample(pool, 5), j] = 1 }   # 5 heterozygous carriers
Gr[, 11] = 0; Gr[sample(N, 2), 11] = 2; Gr[sample(N, 1), 11] = 1          # 2 homozygous + 1 het
Gr = cbind(Gr, Gs[, 1:4])
add("rare_snps_factor3_mixed_df", mk(ind_kernel(ind), ind, cnd3, Gr, lapply(1:3, function(k) rnorm(N) + 0.3 * u)))
add("rare_snps_factor2", mk(ind_kernel(ind), ind, cnd, Gr, lapply(1:3, function(k) rnorm(N) + 0.3 * u)))
# 8 collinear designs: between-individual condition (GRM, 1 sample per individual), SNP whose carriers are all
#   treated (G x C = G), SNP equal to the condition (X1 singular), duplicated SNPs and a sign-flipped copy
ni2 = 200; ind2 = sprintf("j%03d", 1:ni2); cb = rep(c("a", "b"), each = ni2 / 2); Kg2 = grm(ni2, ind2)
Gc = geno(ni2, 10); Gc[, 1] = 0; Gc[sample(which(cb == "b"), 20), 1] = 1; Gc[, 2] = 2 * (cb == "b")
Gc[, 3] = Gc[, 4]; Gc[, 5] = 2 - Gc[, 4]
add("collinear_interaction_and_snps", mk(Kg2, ind2, cb, Gc, lapply(1:3, function(k) 0.5 * (cb == "b") + rnorm(ni2))), "report")
# 9 no condition, no covariates (p0 = 0, q = 0)
add("eqtl_only_p0_0", mk(Kg2, ind2, rep(1, ni2), geno(ni2, 25), lapply(1:3, function(k) 0.4 * Gc[, 6] * (k == 1) + rnorm(ni2))))
# 10 covariates: duplicated covariate (X0 rank-deficient), and a raw-scale (unstandardized) covariate
cv = matrix(rnorm(N * 3), N); cv = cbind(cv, cv[, 1])
add("duplicated_covariate", mk(ind_kernel(ind), ind, cnd, Gs[, 1:10], lapply(1:3, function(k) rnorm(N) + 0.3 * u), covs = cv), "report")
add("raw_scale_covariate", mk(ind_kernel(ind), ind, cnd, Gs[, 1:10], lapply(1:3, function(k) rnorm(N) + 0.3 * u),
                              covs = cbind(age = 1e4 + rnorm(N) * 1e3, rnorm(N) * 1e-4), raw_X0 = TRUE))
# 11 near-perfect fit (y almost in the span of the condition): s -> 0, cancellation in y'Wy - rhs'b
add("near_perfect_fit", mk(ind_kernel(ind), ind, cnd, Gs[, 1:10], lapply(c(1e-3, 1e-6, 1e-9), function(e) (cnd == "b") + e * rnorm(N))), "report")
# 12 SNP counts at block / SIMD boundaries, small and odd n
for (S in c(1, 2, 7, 64, 65, 257)) add(sprintf("n_snps_%d", S), mk(ind_kernel(ind), ind, cnd, geno(N, S), list(rnorm(N) + 0.3 * u)))
ns = 31; inds = sprintf("s%02d", 1:ns)
add("small_odd_n31", mk(grm(ns, inds), inds, rnorm(ns), geno(ns, 9, rep(0.3, 9)), lapply(1:3, function(k) rnorm(ns))))
# 13 boot = FALSE path
P = cases$default_kernel_balanced; P$boot = FALSE; add("no_bootstrap", P)

options(width = 200)
tol = c(l = 1e-6, rel = 1e-6)
ext = NULL
if (startsWith(cand, "lib:")) { fi = tempfile(fileext = ".rds"); fo = tempfile(fileext = ".rds"); saveRDS(cases, fi)
  system2("Rscript", c("bench/proto/review/run_lib.R", sub("^lib:", "", cand), fi, fo)); ext = readRDS(fo) }
res = do.call(rbind, lapply(names(cases), function(nm) {
  P = cases[[nm]]
  run1 = function(f) tryCatch(run_impl(f, P), error = function(e) e)
  b = run1(baseline_fit); c = if (is.null(ext)) run1(fit_gene_proto) else ext[[nm]]
  if (inherits(b, "error") || inherits(c, "error")) {
    same = inherits(b, "error") && inherits(c, "error")
    return(data.frame(case = nm, kind = tag[[nm]], max_dl = NA, rel_beta = NA, rel_se = NA, df = NA, na_mismatch = NA,
                      result = if (same) "PASS (both error)" else paste("FAIL error:", conditionMessage(if (inherits(c, "error")) c else b))))
  }
  x = compare(b, c)
  A = do.call(rbind, lapply(b$out, `[[`, "res")); B = do.call(rbind, lapply(c$out, `[[`, "res"))
  relf = function(j) { d = abs(A[, j] - B[, j]); m = pmax(abs(A[, j]), abs(B[, j]), if (j == 4) A[, 5] else 0);   # beta: relative to max(|beta|, se)
    max(c(0, ifelse(d == 0, 0, d / m)), na.rm = TRUE) }
  dl = max(x[c("l0", "l_geno", "l_interact", "l_boot_geno", "l_boot_interact")], na.rm = TRUE)
  ok = (is.infinite(dl) || dl <= tol[["l"]]) && relf(4) <= tol[["rel"]] && relf(5) <= tol[["rel"]] &&
    identical(A[, 3], B[, 3]) && x[["na_mismatch"]] == 0
  data.frame(case = nm, kind = tag[[nm]], max_dl = dl, rel_beta = relf(4), rel_se = relf(5),
             df = identical(A[, 3], B[, 3]), na_mismatch = x[["na_mismatch"]],
             result = if (ok) "PASS" else if (tag[[nm]] == "report") "DIFF (report)" else "FAIL")
}))
cat("candidate:", cand, "\n"); print(format(res, digits = 3), row.names = FALSE)
cat(sprintf("%d PASS, %d FAIL, %d DIFF on report-only cases\n", sum(startsWith(res$result, "PASS")),
            sum(startsWith(res$result, "FAIL")), sum(startsWith(res$result, "DIFF"))))
