# Kernel-level harness: the per-gene inputs of suez_fit_gene (src/lmm.cpp) for a benchmark dataset, so a
# prototype of the C++ fitter can be timed and compared with the baseline in isolation.
#   source("bench/proto/harness.R")
#   P = prepare("M", n_genes = 50)                     # rotation, design and per-gene genotypes
#   base = run_impl(baseline_fit, P)                   # list of per-gene results + time
#   Rcpp::sourceCpp("bench/proto/p1_gemm_grid.cpp")    # defines fit_gene_proto()
#   new = run_impl(fit_gene_proto, P); compare(base, new)
suppressPackageStartupMessages({ library(suez, lib.loc = "bench/lib/baseline"); library(Rcpp) })

prepare = function(dsn, n_genes = Inf, boot = TRUE) {
  d = readRDS(file.path("bench/data", paste0(dsn, ".rds")))
  input = suez:::quantile_normalize(d$input); anno = d$anno
  K = d$sample_kernel
  if (is.null(K)) { K = outer(anno$individual, anno$individual, "=="); class(K) = "numeric" }
  eig = suez:::kernel_contrast_eigen(K)
  one = length(unique(anno$condition)) == 1
  C = if (one) matrix(0, nrow(anno), 0) else
    if (is.numeric(anno$condition)) cbind(anno$condition) else model.matrix(~ condition, data = anno)[, -1, drop = FALSE]
  X0 = if (one) matrix(0, nrow(anno), 0) else model.matrix(~ condition, data = anno)
  if (!is.null(d$covariates)) X0 = cbind(X0, d$covariates)
  X0 = suez:::standardize_cols(X0[, apply(X0, 2, sd) > 0, drop = FALSE])
  genes = head(d$geneloc$geneid, n_genes)
  set.seed(42); seeds = sample.int(.Machine$integer.max, length(genes))
  per = lapply(seq_along(genes), function(i) {
    gl = d$geneloc[d$geneloc$geneid == genes[i], ]
    sel = d$snploc$chr == gl$chr & d$snploc$pos > gl$left - 1e5 & d$snploc$pos < gl$right + 1e5
    g = d$genotype[as.character(d$snploc$snpid[sel]), , drop = FALSE]
    g = if (anyNA(g)) suez:::easy_impute(g) else { storage.mode(g) = "integer"; g }
    G = t(g[, anno$individual, drop = FALSE]); storage.mode(G) = "double"
    G = G[, rowSums(g) >= 5 & apply(G, 2, sd) > 0, drop = FALSE]
    list(gene = genes[i], seed = seeds[i], Uy = drop(crossprod(eig$vectors, input[genes[i], ])), G = G)
  })
  list(U = eig$vectors, lambda = eig$values, UX0 = crossprod(eig$vectors, X0), C = C, boot = boot, genes = per)
}

baseline_fit = function(U, lambda, Uy, UX0, G, C, boot) .Call("suez_fit_gene", U, lambda, Uy, UX0, G, C, boot, PACKAGE = "suez")

run_impl = function(fun, P, reps = 1) {
  times = numeric(reps)
  for (r in seq_len(reps)) {
    t0 = proc.time()[["elapsed"]]
    out = lapply(P$genes, function(g) { set.seed(g$seed); fun(P$U, P$lambda, g$Uy, P$UX0, g$G, P$C, P$boot) })
    times[r] = proc.time()[["elapsed"]] - t0
  }
  list(out = out, time = times, n_snp = sum(sapply(P$genes, function(g) ncol(g$G))))
}

# max |difference| per output column (l_geno, l_interact, df, beta, se, l_boot_geno, l_boot_interact)
# and the number of NA mismatches
compare = function(a, b) {
  A = do.call(rbind, lapply(a$out, `[[`, "res")); B = do.call(rbind, lapply(b$out, `[[`, "res"))
  l0 = max(abs(sapply(a$out, `[[`, "l0") - sapply(b$out, `[[`, "l0")))
  d = sapply(seq_len(ncol(A)), function(j) suppressWarnings(max(abs(A[, j] - B[, j]), na.rm = TRUE)))
  rel = sapply(4:5, function(j) suppressWarnings(max(abs(A[, j] - B[, j]) / abs(A[, j]), na.rm = TRUE)))
  c(l0 = l0, setNames(d, c("l_geno", "l_interact", "df", "beta", "se", "l_boot_geno", "l_boot_interact")),
    rel_beta = rel[1], rel_se = rel[2], na_mismatch = sum(is.na(A) != is.na(B)))
}
