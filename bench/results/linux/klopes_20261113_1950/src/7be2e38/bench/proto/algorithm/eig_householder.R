# Idea 3: kernel_contrast_eigen with the Householder reflection applied in O(n^2) instead of dense O(n^3) products.
suppressPackageStartupMessages(library(suez, lib.loc = "bench/lib/baseline"))
kce_fast = function(K) {
  K = as.matrix(K); n = nrow(K)
  stopifnot(nrow(K) == ncol(K), isSymmetric(unname(K), tol = 1e-8))
  H1 = qr.Q(qr(matrix(1, n, 1)))[, 1]           # first column of the complete Q = H e1, H = I - 2 u u'
  u = -H1; u[1] = u[1] + 1; u = u / sqrt(sum(u^2))
  Ku = drop(K %*% u); uKu = sum(u * Ku)
  KH = K - 2 * outer(Ku, u)                        # K H
  HKH = KH - 2 * outer(u, drop(crossprod(u, KH)))  # H K H
  e = eigen(HKH[-1, -1], symmetric = TRUE)
  if (min(e$values) < -1e-6 * max(abs(e$values))) stop("not PSD")
  V = rbind(0, e$vectors)
  list(vectors = V - 2 * outer(u, drop(crossprod(u, V))), values = pmax(e$values, 0))
}
for (dsn in c("S", "M_grm", "L", "XL")) {
  d = readRDS(file.path("bench/data", paste0(dsn, ".rds"))); anno = d$anno
  K = d$sample_kernel
  if (is.null(K)) { K = outer(anno$individual, anno$individual, "=="); class(K) = "numeric" }
  n = nrow(K)
  t0 = proc.time()[[3]]; Q = qr.Q(qr(matrix(1, n, 1)), complete = TRUE)[, -1, drop = FALSE]
  t1 = proc.time()[[3]]; M = crossprod(Q, K %*% Q)
  t2 = proc.time()[[3]]; e = eigen(M, symmetric = TRUE)
  t3 = proc.time()[[3]]; V = Q %*% e$vectors
  t4 = proc.time()[[3]]; f = kce_fast(K); t5 = proc.time()[[3]]
  b = suez:::kernel_contrast_eigen(K)
  Mf = { H1 = Q } ; rm(Mf)
  cat(sprintf("%s n=%d  base: qr.Q %.2fs  Q'KQ %.2fs  eigen %.2fs  QV %.2fs | fast total %.2fs | max|dlambda| %.2g  max|dU| %.2g\n",
      dsn, n, t1 - t0, t2 - t1, t3 - t2, t4 - t3, t5 - t4, max(abs(f$values - b$values)), max(abs(f$vectors - b$vectors))))
}
