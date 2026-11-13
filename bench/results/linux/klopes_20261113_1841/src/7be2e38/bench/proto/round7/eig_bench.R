# Round 7: cost of kernel_contrast_eigen variants at N = 1000, 2000 (GRM-like and default kernels).
suppressPackageStartupMessages(library(suez, lib.loc = "bench/lib/r6_p9"))
Rcpp::sourceCpp("bench/proto/round7/eig_solvers.cpp")
source("bench/proto/algorithm/eig_householder_fun.R")
hh_proj = function(K) {   # Householder projection only (O(n^2)): H K H without the first row/column, and u
  n = nrow(K); H1 = qr.Q(qr(matrix(1, n, 1)))[, 1]
  u = -H1; u[1] = u[1] + 1; u = u / sqrt(sum(u^2))
  KH = K - 2 * outer(drop(K %*% u), u); HKH = KH - 2 * outer(u, drop(crossprod(u, KH)))
  list(A = HKH[-1, -1], u = u) }
hh_back = function(u, V) { V = rbind(0, V); V - 2 * outer(u, drop(crossprod(u, V))) }
for (n in c(1000, 2000)) {
  set.seed(1); Z = matrix(rnorm(n * 300), n); K = tcrossprod(Z) / 300
  t_base = system.time(e0 <- suez:::kernel_contrast_eigen(K))[["elapsed"]]
  t_proj = system.time(p <- hh_proj(K))[["elapsed"]]
  t_R = system.time(eR <- eigen(p$A, symmetric = TRUE))[["elapsed"]]
  t_E = system.time(eE <- eigen_sym(p$A))[["elapsed"]]
  t_back = system.time(V <- hh_back(p$u, eE$vectors))[["elapsed"]]
  # same subspaces? compare projectors on the top 10 (distinct) eigenvalues and values
  cat(sprintf("n=%d base %.2fs | householder proj %.2fs, R eigen %.2fs, Eigen solver %.2fs, back %.2fs | max dval %.2g, top-10 projector diff %.2g\n",
      n, t_base, t_proj, t_R, t_E, t_back, max(abs(e0$values - eE$values)),
      max(abs(tcrossprod(e0$vectors[, 1:10]) - tcrossprod(V[, 1:10])))))
}
