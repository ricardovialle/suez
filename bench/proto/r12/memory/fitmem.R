# Rscript fitmem.R <lib> <dataset> <n_genes>: footprint (now / peak so far) around set-up and each suez_fit_gene call
a = commandArgs(TRUE); suppressPackageStartupMessages(library(suez, lib.loc = a[1])); ns = asNamespace("suez")
fp = function(s) { o = system(sprintf("footprint -p %d --noCategories -f bytes", Sys.getpid()), intern = TRUE)
  v = as.numeric(sub(".*: ([0-9]+) B.*", "\\1", grep("phys_footprint", o, value = TRUE))) / 2^20
  cat(sprintf("%-30s footprint %7.1f MB   peak %7.1f MB\n", s, v[1], v[2])) }
d = readRDS(file.path("/Users/rvialle/suez/bench/data", paste0(a[2], ".rds"))); anno = d$anno
K = d$sample_kernel; if (is.null(K)) { Z = matrix(0, nrow(anno), length(unique(anno$individual))); Z[cbind(seq_len(nrow(anno)), match(anno$individual, unique(anno$individual)))] = 1; K = tcrossprod(Z); rm(Z) }
eig = ns$kernel_contrast_eigen(K); rm(K); U = eig$vectors; lambda = eig$values; invisible(gc()); fp("set-up done")
C = model.matrix(~ condition, data = anno)[, -1, drop = FALSE]; X0 = ns$standardize_cols(C); UX0 = crossprod(U, X0)
input = ns$quantile_normalize(d$input); snp = d$snploc; gl = d$geneloc
for (j in seq_len(as.integer(a[3]))) { gl1 = gl[j, ]; s = as.character(snp$snpid[snp$chr == gl1$chr & snp$pos > gl1$left - 1e5 & snp$pos < gl1$right + 1e5])
  G = t(round(d$genotype[s, anno$individual, drop = FALSE])); G = G[, colSums(G) >= 5, drop = FALSE]
  Uy = drop(crossprod(U, input[j, ])); if (nzchar(Sys.getenv("GCG"))) cat("gc s", system.time(gc(FALSE))[[3]], "\n"); fp(sprintf("gene %d before (%d SNPs)", j, ncol(G)))
  f = .Call(ns$suez_fit_gene, U, lambda, Uy, UX0, G, C, TRUE); rm(G); fp(sprintf("gene %d after", j)) }
