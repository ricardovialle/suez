# Eigen-solver timings for kernel_contrast_eigen-sized problems (GRM-like K, N-1 x N-1): Eigen (base12), R eigen()
# (reference LAPACK dsyevr here), Accelerate dsyevd / dsyevr. Usage: Rscript eig_bench.R N ...
Sys.setenv(PKG_LIBS = "-framework Accelerate"); Rcpp::sourceCpp("bench/proto/r12/numerical/eig/eig_acc.cpp")
for (n in as.integer(commandArgs(TRUE))) {
  set.seed(1); Z = scale(matrix(rbinom(n * 3000, 2, 0.3), n)); K = tcrossprod(Z) / 3000; rm(Z)
  K = K[-1, -1]
  tm = c(); ev = list()
  for (m in c(1, 2, 0)) { t = system.time(ev[[m + 1]] <- eig_acc(K, m)); tm[c("eigen", "acc_dsyevd", "acc_dsyevr")[m + 1]] = t[["elapsed"]] }
  if (n <= 2000) { t = system.time(e = eigen(K, symmetric = TRUE)); tm["R_eigen_ref"] = t[["elapsed"]] }
  v0 = ev[[1]]$values
  dv = sapply(ev[2:3], function(e) max(abs(e$values - v0)) / max(v0))
  # projector on the top 10 eigenvectors (distinct eigenvalues): basis-invariant check
  top = function(e) { i = order(e$values, decreasing = TRUE)[1:10]; tcrossprod(e$vectors[, i]) }
  dp = sapply(ev[2:3], function(e) max(abs(top(e) - top(ev[[1]]))))
  res = sapply(ev, function(e) { V = e$vectors[, 1:50]; max(abs(K %*% V - V %*% diag(e$values[1:50]))) })
  cat(sprintf("N=%d  %s  | max rel dvalues vs Eigen %s  top-10 projector diff %s  residual(50 vecs) %s\n", n + 1,
      paste(names(tm), sprintf("%.2fs", tm), collapse = " "), paste(signif(dv, 2), collapse = "/"), paste(signif(dp, 2), collapse = "/"), paste(signif(res, 2), collapse = "/")))
}
