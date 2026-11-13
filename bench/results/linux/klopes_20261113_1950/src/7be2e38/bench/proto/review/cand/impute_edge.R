# easy_impute edge cases with one installed build: Rscript impute_edge.R <lib dir>
a = commandArgs(TRUE); suppressPackageStartupMessages(library(suez, lib.loc = a[1]))
set.seed(3); ni = 60; G = matrix(rbinom(8 * ni, 2, 0.3), 8, ni, dimnames = list(paste0("rs", 1:8), paste0("I", 1:ni)))
G[cbind(sample(8, 20, TRUE), sample(ni, 20, TRUE))] = NA
try1 = function(nm, g) { r = tryCatch(suez:::easy_impute(g), error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-38s %s\n", nm, if (is.character(r)) r else sprintf("ok, NA left %d, sum %d", sum(is.na(r)), sum(r, na.rm = TRUE)))) }
try1("typical (8 SNPs, 20 NA)", G)
g = G; g[3, ] = 0L; g[3, 1:3] = NA; try1("monomorphic SNP with missing calls", g)
g = G; g[4, ] = NA; try1("SNP with all calls missing", g)
g = G; g[5, ] = NA; g[5, 7] = 1L; try1("SNP with one non-missing call", g)
try1("one SNP", G[1, , drop = FALSE])
g = G[rep(1, 8), ]; g[2, 5] = NA; try1("8 identical SNPs (k = 0)", g)
try1("more SNPs than individuals", matrix(c(rbinom(80 * 30, 2, 0.3)), 80, 30) |> (\(m) { m[1:5] = NA; m })())
# end to end: a gene whose cis window holds such a SNP (debug = FALSE drops the failing gene silently)
Y = matrix(rnorm(2 * ni), 2, dimnames = list(c("gA", "gB"), paste0("s", 1:ni)))
g = G; g[3, ] = 0L; g[3, 1:3] = NA
r = map_interaction_qtl(Y, g, data.frame(geneid = c("gA", "gB"), chr = "1", left = c(100, 1e6), right = c(110, 1e6 + 10)),
  data.frame(snpid = rownames(g), chr = "1", pos = c(101:104, 1e6 + 1:4)), data.frame(individual = paste0("I", 1:ni), condition = 1),
  permutation_approach = "none", cisdist = 100)
cat("end to end: genes in result:", paste(unique(r$gene), collapse = ","), "(gA holds the monomorphic SNP)\n")
# proposed fix: zero every NaN of the scaled matrix (as 0.4: constant rows give 0/0), keep miss for the fill-in
fixed = function(geno, prop_var = 0.95) {
  temp = t(scale(t(geno))); miss = is.na(geno); temp[is.na(temp)] = 0
  wide = nrow(temp) <= ncol(temp); e = .Call(suez:::suez_eigen_sym, if (wide) tcrossprod(temp) else crossprod(temp))
  v = pmax(e$values, 0); v = v / sum(v); W = e$vectors[, seq_len(sum(cumsum(v) < prop_var)), drop = FALSE]
  recon = if (wide) W %*% crossprod(W, temp) else (temp %*% W) %*% t(W)
  temp[miss] = recon[miss]; temp = suez:::unscale(temp); temp = round(temp); class(temp) = "integer"; temp
}
g2 = G; g2[5, ] = NA; g2[5, 7] = 1L
for (x in list(G, g, g2)) cat("fixed version: NA left", sum(is.na(fixed(x))), " sum", sum(fixed(x), na.rm = TRUE), "\n")
