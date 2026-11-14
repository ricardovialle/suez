source("bench/proto/r12/algorithm/zz_basis.R")
chk = function(ind) {
  b = zz_contrast_basis(ind); U = sp_dense(b$U); N = length(ind)
  Z = outer(ind, unique(ind), "==") * 1; K = tcrossprod(Z)
  c(N = N, orth = max(abs(crossprod(U) - diag(N - 1))), one = max(abs(colSums(U))),
    eig = max(abs(crossprod(U, K %*% U) - diag(b$values, N - 1))), nnz = length(b$U$x))
}
set.seed(1)
print(chk(rep(1:50, each = 2)))
print(chk(sample(rep(1:40, times = sample(1:4, 40, TRUE)))))
print(chk(c(1:5, 1:5, 6:7)))
d = readRDS("bench/data/M_unbal.rds"); print(chk(d$anno$individual))
for (n in c(1000, 2000)) { ind = rep(seq_len(n), each = 2); print(system.time(b <- zz_contrast_basis(ind))) }
