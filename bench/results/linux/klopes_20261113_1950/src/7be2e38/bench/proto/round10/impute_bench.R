# Round 10: easy_impute via the eigen-decomposition of the smaller Gram matrix vs svd(), on every gene of a
# dataset with missing genotypes. Output must be identical (imputed values are rounded to integers).
suppressPackageStartupMessages(library(suez, lib.loc = "bench/lib/r8_fix"))
eig = function(A) .Call("suez_eigen_sym", A, PACKAGE = "suez")
easy_impute2 = function(geno, prop_var = 0.95) {
  temp = t(scale(t(geno))); miss = is.na(geno); temp[miss] = 0
  # temp = U D V': the rank-k reconstruction U_k D_k V_k' is the projection of temp on U_k (rows <= cols) or
  # of temp' on V_k, so only the eigen-decomposition of the smaller Gram matrix is needed
  wide = nrow(temp) <= ncol(temp)
  e = eig(if (wide) tcrossprod(temp) else crossprod(temp))
  v = pmax(e$values, 0); v = v / sum(v)
  k = sum(cumsum(v) < prop_var)
  W = e$vectors[, seq_len(k), drop = FALSE]
  recon = if (wide) W %*% crossprod(W, temp) else (temp %*% W) %*% t(W)
  temp[miss] = recon[miss]
  temp = suez:::unscale(temp)
  stopifnot(max(abs(temp[!miss] - geno[!miss])) < 1e-10)
  temp = round(temp); class(temp) = "integer"; temp
}
d = readRDS("bench/data/M_miss.rds")
genes = d$geneloc$geneid
win = lapply(genes, function(g) { gl = d$geneloc[d$geneloc$geneid == g, ]
  sel = d$snploc$chr == gl$chr & d$snploc$pos > gl$left - 1e5 & d$snploc$pos < gl$right + 1e5
  d$genotype[as.character(d$snploc$snpid[sel]), , drop = FALSE] })
t1 = system.time(a <- lapply(win, suez:::easy_impute))[["elapsed"]]
t2 = system.time(b <- lapply(win, easy_impute2))[["elapsed"]]
nd = sum(mapply(function(x, y) sum(x != y), a, b))
cat(sprintf("easy_impute %.2fs, Gram-eigen %.2fs (%.1fx), %d of %d genes, imputed entries differing: %d of %d\n",
  t1, t2, t1 / t2, length(win), length(win), nd, sum(sapply(win, function(x) sum(is.na(x))))))
# wide/tall shapes and a small window
set.seed(3); for (dims in list(c(30, 200), c(400, 200), c(5, 50))) {
  g = matrix(rbinom(prod(dims), 2, 0.3), dims[1]); g[sample(length(g), length(g) / 50)] = NA
  cat(dims, "identical:", identical(suez:::easy_impute(g), easy_impute2(g)), "\n") }
