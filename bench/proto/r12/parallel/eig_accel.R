# Eigen's solver vs Accelerate dsyevd/dsyevr on the Householder-projected kernel (as kernel_contrast_eigen).
#   VECLIB_MAXIMUM_THREADS=1 Rscript bench/proto/r12/parallel/eig_accel.R <threads label> [N ...]
a = commandArgs(TRUE); lab = a[1]; Ns = if (length(a) > 1) as.integer(a[-1]) else c(1000, 2000, 4000)
Sys.setenv(PKG_LIBS = "-framework Accelerate")
Rcpp::sourceCpp("bench/proto/r12/parallel/eig_accel.cpp", cacheDir = file.path(tempdir(), "eigacc"))
proj = function(K) { n = nrow(K); u = -qr.Q(qr(matrix(1, n, 1)))[, 1]; u[1] = u[1] + 1; u = u / sqrt(sum(u^2))
  K = K - 2 * outer(drop(K %*% u), u); K = K - 2 * outer(u, drop(crossprod(u, K))); K[-1, -1] }
set.seed(1)
for (N in Ns) for (kind in c("grm", "default")) {
  if (kind == "grm") { X = matrix(rbinom(N * 3000, 2, 0.3), N); X = scale(X); K = tcrossprod(X) / ncol(X); rm(X) }
  else { ind = rep(seq_len(N / 2), each = 2); K = outer(ind, ind, "=="); class(K) = "numeric" }
  A = proj(K); rm(K)
  te = if (N <= 2000) system.time(e1 <- eig_eigen(A))[["elapsed"]] else NA
  td = system.time(e2 <- eig_accel(A, FALSE))[["elapsed"]]
  tr = system.time(e3 <- eig_accel(A, TRUE))[["elapsed"]]
  dv = if (!is.na(te)) max(abs(e1$values - e2$values)) else NA
  # projector on the leading (distinct for grm) eigenvector
  dp = if (!is.na(te) && kind == "grm") max(abs(tcrossprod(e1$vectors[, 1]) - tcrossprod(e2$vectors[, 1]))) else NA
  cat(sprintf("%s\t%s\t%d\teigen=%.2f\tdsyevd=%.2f\tdsyevr=%.2f\tdval=%.1e\tdproj1=%.1e\n", lab, kind, N, te, td, tr, dv, dp))
}
