# kernel_contrast_eigen edge cases with one installed build: Rscript eig_edge.R <lib dir>
a = commandArgs(TRUE); suppressPackageStartupMessages(library(suez, lib.loc = a[1])); f = suez:::kernel_contrast_eigen
t1 = function(nm, expr) { r = tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e))); cat(sprintf("%-34s %s\n", nm, if (is.character(r)) r else r)) }
chk = function(K) { e = f(K); n = nrow(K); U = e$vectors
  sprintf("ok: n-1 = %d values, decreasing %s, |U'1| %.1e, |U'U - I| %.1e, |U diag U' - PKP| %.1e", length(e$values), !is.unsorted(rev(e$values)),
    max(abs(crossprod(U, rep(1, n)))), max(abs(crossprod(U) - diag(ncol(U)))),
    { P = diag(n) - 1 / n; max(abs(U %*% (e$values * t(U)) - P %*% K %*% P)) }) }
set.seed(2); A = matrix(rnorm(400), 20); K = tcrossprod(A) / 20
t1("random PSD 20 x 20", chk(K))
t1("n = 2", chk(matrix(c(2, 1, 1, 2), 2)))
t1("n = 3, default kernel", chk(outer(c(1, 1, 2), c(1, 1, 2), "==") * 1))
t1("integer storage kernel", chk(outer(1:6 %% 3, 1:6 %% 3, "==") + 0L))
t1("not PSD", chk(K - 3 * diag(20)))
t1("asymmetric", chk(K + outer(1:20, rep(1e-3, 20))))
K2 = K; K2[3, 4] = K2[4, 3] = NA; t1("NA entry", chk(K2))
K3 = K; K3[3, 4] = K3[4, 3] = Inf; t1("Inf entry", chk(K3))
cat("values vs eigen(): max diff", max(abs(f(K)$values - { Q = qr.Q(qr(matrix(1, 20, 1)), complete = TRUE)[, -1]; eigen(crossprod(Q, K %*% Q), symmetric = TRUE)$values })), "\n")
