# fixed vs per-SNP cost of suez_fit_gene at N = 200 (P_many shapes)
suppressPackageStartupMessages(library(suez, lib.loc = "bench/lib/base12")); ns = asNamespace("suez")
d = readRDS("bench/data/P_many.rds"); anno = d$anno
eig = ns$kernel_contrast_eigen(tcrossprod(outer(anno$individual, unique(anno$individual), "==") * 1)); U = eig$vectors; lambda = eig$values
C = model.matrix(~ condition, data = anno)[, -1, drop = FALSE]; X0 = ns$standardize_cols(C); UX0 = crossprod(U, X0)
set.seed(1); Gall = t(d$genotype[sample(nrow(d$genotype), 200), anno$individual]); storage.mode(Gall) = "double"
Gall = Gall[, colSums(Gall) >= 5 & apply(Gall, 2, sd) > 0]
y = crossprod(U, rnorm(nrow(anno)))
for (boot in c(FALSE, TRUE)) for (S in c(1, 2, 10, 50)) {
  G = Gall[, 1:S, drop = FALSE]; r = 200 / S
  t = system.time(for (k in 1:r) .Call(ns$suez_fit_gene, U, lambda, drop(y), UX0, G, C, boot))[[3]]
  cat(sprintf("boot %d S %3d: %.3f ms per call, %.3f ms per SNP\n", boot, S, 1e3 * t / r, 1e3 * t / 200)) }
# real P_many genes
input = ns$quantile_normalize(d$input[1:200, ]); snp = d$snploc; gl = d$geneloc
tt = 0; ns_ = 0
for (j in 1:200) { gl1 = gl[j, ]; s = as.character(snp$snpid[snp$chr == gl1$chr & snp$pos > gl1$left - 1e5 & snp$pos < gl1$right + 1e5])
  G = t(d$genotype[s, anno$individual, drop = FALSE]); storage.mode(G) = "double"; G = G[, colSums(G) >= 5 & apply(G, 2, sd) > 0, drop = FALSE]
  Uy = drop(crossprod(U, input[j, ])); ns_ = ns_ + ncol(G)
  tt = tt + system.time(.Call(ns$suez_fit_gene, U, lambda, Uy, UX0, G, C, TRUE))[[3]] }
cat(sprintf("real genes: %.3f ms per gene, %.1f SNPs per gene\n", 1e3 * tt / 200, ns_ / 200))
y2 = drop(crossprod(U, input[1, ])); G = Gall[, 1:10]
cat(sprintf("random G, real y: %.3f ms\n", 1e3 * system.time(for (k in 1:50) .Call(ns$suez_fit_gene, U, lambda, y2, UX0, G, C, TRUE))[[3]] / 50))
