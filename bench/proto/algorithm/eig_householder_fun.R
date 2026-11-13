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
