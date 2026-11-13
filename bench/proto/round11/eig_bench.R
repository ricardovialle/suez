# Round 11: kernel_contrast_eigen under the BLAS/LAPACK this R links to.
# Usage: Rscript eig_bench.R <suez lib> <out.tsv> [N ...]
# Compares, on a GRM-like kernel (Z Z' / m, m = 3000 SNPs):
#   eigen_cpp : suez's kernel_contrast_eigen (Householder + Eigen's SelfAdjointEigenSolver, suez >= 0.5)
#   lapack    : the same Householder projection + R's eigen() (LAPACK dsyevr through R's BLAS/LAPACK)
#   dense     : suez 0.4's Q'KQ / QV dense products + eigen() (reference point)
# and checks the eigenvalues and the projector onto the leading eigenvectors agree.
args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]; out <- args[2]
Ns <- if (length(args) > 2) as.integer(args[-(1:2)]) else c(500L, 1000L, 2000L, 5000L)
library(suez, lib.loc = lib)
kce <- utils::getFromNamespace("kernel_contrast_eigen", "suez")

householder_lapack <- function(K) {          # kernel_contrast_eigen with eigen() in place of suez_eigen_sym
  n <- nrow(K)
  u <- -qr.Q(qr(matrix(1, n, 1)))[, 1]; u[1] <- u[1] + 1; u <- u / sqrt(sum(u^2))
  K <- K - 2 * outer(drop(K %*% u), u)
  K <- K - 2 * outer(u, drop(crossprod(u, K)))
  e <- eigen(K[-1, -1, drop = FALSE], symmetric = TRUE)
  V <- rbind(0, e$vectors)
  list(vectors = V - 2 * outer(u, drop(crossprod(u, V))), values = pmax(e$values, 0))
}
dense <- function(K) {
  n <- nrow(K); Q <- qr.Q(qr(matrix(1, n, 1)), complete = TRUE)[, -1, drop = FALSE]
  e <- eigen(crossprod(Q, K %*% Q), symmetric = TRUE)
  list(vectors = Q %*% e$vectors, values = pmax(e$values, 0))
}

blas <- tryCatch(extSoftVersion()[["BLAS"]], error = function(e) NA)
info <- c(R = R.version.string, BLAS = blas, LAPACK = La_library(), suez = as.character(packageVersion("suez", lib.loc = lib)),
          threads = Sys.getenv("OPENBLAS_NUM_THREADS", "unset"), sysname = Sys.info()[["sysname"]], machine = Sys.info()[["machine"]])
print(info)
res <- list()
for (n in Ns) {
  set.seed(n); Z <- scale(matrix(rbinom(n * 3000, 2, 0.3), n)); K <- tcrossprod(Z) / 3000; rm(Z)
  t <- c(eigen_cpp = system.time(a <- kce(K))[["elapsed"]],
         lapack = system.time(b <- householder_lapack(K))[["elapsed"]],
         dense = if (n <= 2000) system.time(d <- dense(K))[["elapsed"]] else NA)
  k <- 20   # leading eigenvalues are distinct for a GRM: compare their projector
  P <- function(e) tcrossprod(e$vectors[, 1:k])
  res[[length(res) + 1]] <- data.frame(N = n, t(t), dvalues = max(abs(a$values - b$values)) / max(a$values),
                                       dproj = max(abs(P(a) - P(b))), as.list(info))
  print(res[[length(res)]][, 1:6])
  rm(K, a, b); gc()
}
write.table(do.call(rbind, res), out, sep = "\t", quote = FALSE, row.names = FALSE)
