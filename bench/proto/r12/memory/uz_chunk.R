Rcpp::sourceCpp("uz_chunk.cpp")
set.seed(1); N = 2000; U = qr.Q(qr(matrix(rnorm(N * N), N)))[, -1]
G = matrix(rbinom(N / 2 * 256, 2, 0.3), N / 2)[rep(1:(N / 2), each = 2), ]; storage.mode(G) = "double"
for (rb in c(256, 512, 1024)) { r = uz(U, G, rb); cat("rows per block", rb, "nc", r[[3]], "identical", identical(r[[1]], r[[2]]), "max diff", max(abs(r[[1]] - r[[2]])), "\n") }
