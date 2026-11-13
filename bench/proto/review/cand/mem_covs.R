# Peak RAM of one fit_gene call with many covariates, a GRM (many eigenvalue groups) and a large cis window.
#   /usr/bin/time -l Rscript bench/proto/review/cand/mem_covs.R <impl.cpp|baseline> <n_cov> <S>
source("bench/proto/harness.R"); a = commandArgs(TRUE)
f = if (a[1] == "baseline") baseline_fit else { e = new.env(); Rcpp::sourceCpp(a[1], env = e); e$fit_gene_proto }
nc = as.integer(a[2]); S = as.integer(a[3]); set.seed(1); ni = 200; ind = rep(1:ni, each = 2); N = 2 * ni
Z = matrix(rbinom(ni * 2000, 2, 0.3), ni); Z = scale(Z[, apply(Z, 2, sd) > 0]); K = (tcrossprod(Z) / ncol(Z))[ind, ind]
eig = suez:::kernel_contrast_eigen(K); cond = rep(c(0, 1), ni)
X0 = suez:::standardize_cols(cbind(cond, matrix(rnorm(N * nc), N))); G = sapply(runif(S, 0.1, 0.5), function(p) rbinom(ni, 2, p))[ind, ]
gc(reset = TRUE); m0 = sum(gc()[, 2])
t = system.time(r <- f(eig$vectors, eig$values, drop(crossprod(eig$vectors, rnorm(N))), crossprod(eig$vectors, X0), G, cbind(cond), TRUE))
cat(sprintf("RESULT %s n_cov=%d S=%d elapsed=%.2f\n", basename(a[1]), nc, S, t[["elapsed"]]))
