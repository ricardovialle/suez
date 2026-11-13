Rcpp::sourceCpp("gemm_colind.cpp")
set.seed(1); bad = 0
for (n in c(199, 500, 999, 2000)) for (k in c(100, 500, 1000, 2000)) {
  A = matrix(rnorm(k * n), k); B = matrix(rbinom(k * 700, 2, .3) / 0.7, k)
  full = prodE(A, B)
  for (sub in list(1, 1:3, 5:20, 37:300, 1:700, 650:700)) { if (!identical(full[, sub, drop = FALSE], prodE(A, B[, sub, drop = FALSE]))) { bad = bad + 1; cat("diff n", n, "k", k, "cols", range(sub), max(abs(full[, sub] - prodE(A, B[, sub, drop = FALSE]))), "\n") } }
}
cat("non-identical:", bad, "\n")
