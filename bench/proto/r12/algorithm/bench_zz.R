# Interleaved: base12 (dense Householder+Eigen basis, base engine) vs structured sparse basis + sparse rotation.
# Rscript bench/proto/r12/algorithm/bench_zz.R <dataset> <n_genes> <reps>
a = commandArgs(TRUE); dsn = a[1]; ng = as.numeric(a[2]); reps = as.integer(a[3])
Sys.setenv(SUEZ_HARNESS_LIB = "bench/lib/base12"); source("bench/proto/harness.R")
source("bench/proto/r12/algorithm/zz_basis.R")
Rcpp::sourceCpp(if (length(a) > 3) a[4] else "bench/proto/r12/algorithm/p12_zz.cpp", cacheDir = "bench/proto/r12/algorithm/.cache")
d = readRDS(file.path("bench/data", paste0(dsn, ".rds"))); anno = d$anno
stopifnot(is.null(d$sample_kernel))
# set-up: base kernel + eigen vs analytic basis (interleaved)
ts = sapply(1:reps, function(r) {
  t1 = system.time({ K = outer(anno$individual, anno$individual, "=="); class(K) = "numeric"; e = suez:::kernel_contrast_eigen(K); rm(K, e); gc() })[["elapsed"]]
  t2 = system.time(b <- zz_contrast_basis(anno$individual))[["elapsed"]]
  c(t1, t2) })
P = prepare(dsn, n_genes = ng)
b = zz_contrast_basis(anno$individual)
input = suez:::quantile_normalize(d$input)
one = length(unique(anno$condition)) == 1
X0 = if (one) matrix(0, nrow(anno), 0) else model.matrix(~ condition, data = anno)
if (!is.null(d$covariates)) X0 = cbind(X0, d$covariates)
X0 = suez:::standardize_cols(X0[, apply(X0, 2, sd) > 0, drop = FALSE])
Pz = P; Pz$U = b$U; Pz$lambda = b$values; Pz$UX0 = sp_crossprod(b$U, X0)
for (g in seq_along(Pz$genes)) Pz$genes[[g]]$Uy = drop(sp_crossprod(b$U, cbind(input[P$genes[[g]]$gene, ])))
tf = matrix(NA, reps, 2)
for (r in 1:reps) { x = run_impl(baseline_fit, P); tf[r, 1] = x$time; y = run_impl(fit_gene_sparse, Pz); tf[r, 2] = y$time }
# same zz basis through the dense base engine: isolates the sparse rotation (bootstrap draws comparable)
Pd = Pz; Pd$U = sp_dense(b$U); zd = run_impl(baseline_fit, Pd)
cat(sprintf("%s N=%d genes=%d snps=%d\n", dsn, nrow(anno), length(P$genes), x$n_snp))
cat(sprintf("setup  base %.3f s  zz %.4f s  ratio %.1f\n", median(ts[1, ]), median(ts[2, ]), median(ts[1, ]) / median(ts[2, ])))
cat(sprintf("fit    base %.3f s  zz %.3f s  speed-up %.2fx (median of %d, ratios %s)\n", median(tf[, 1]), median(tf[, 2]),
            median(tf[, 1] / tf[, 2]), reps, paste(round(tf[, 1] / tf[, 2], 2), collapse = " ")))
cat("zz sparse vs base12:\n"); print(signif(compare(x, y), 3))
cat("zz sparse vs zz-basis dense base engine:\n"); print(signif(compare(zd, y), 3))
lr = function(o) { A = do.call(rbind, lapply(o$out, `[[`, "res")); 2 * (A[, 7] - A[, 6]) }
la = lr(x); lb = lr(y); cat("bootstrap LR_int KS p (base vs zz):", signif(suppressWarnings(ks.test(la, lb))$p.value, 3),
  " quantiles base:", signif(quantile(la, c(.5, .9, .99), na.rm = TRUE), 3), " zz:", signif(quantile(lb, c(.5, .9, .99), na.rm = TRUE), 3), "\n")
saveRDS(list(la = la, lb = lb, tf = tf, ts = ts), sprintf("bench/proto/r12/algorithm/zz_%s%s.rds", dsn, if (length(a) > 3) "_combo" else ""))
