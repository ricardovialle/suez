# Round 11: suez_fit_gene time at large N for one installed suez build (run once per build, same seed).
# Usage: Rscript fit_bench.R <suez lib> <out.rds> <N> [genes = 20] [snps = 50] [interaction = 0|1]
# Inputs mimic eqtl_smk's SUEZ engine: GRM-like kernel, intercept + 5 covariates, one sample per donor.
# The kernel eigen-decomposition is done with R's eigen() (identical for every build) and not timed.
args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]; out <- args[2]; n <- as.integer(args[3])
genes <- if (length(args) > 3) as.integer(args[4]) else 20L
snps <- if (length(args) > 4) as.integer(args[5]) else 50L
inter <- length(args) > 5 && args[6] == "1"
library(suez, lib.loc = lib)
sym <- utils::getFromNamespace("suez_fit_gene", "suez")   # native symbol
fit <- function(...) .Call(sym, ...)

set.seed(1)
Z <- scale(matrix(rbinom(n * 3000, 2, 0.3), n)); K <- tcrossprod(Z) / 3000; rm(Z)
Q <- qr.Q(qr(matrix(1, n, 1)), complete = TRUE)[, -1, drop = FALSE]
e <- eigen(crossprod(Q, K %*% Q), symmetric = TRUE); U <- Q %*% e$vectors; lambda <- pmax(e$values, 0); rm(Q, K, e)
X0 <- scale(matrix(rnorm(n * 5), n)); UX0 <- crossprod(U, X0)
C <- if (inter) matrix(rbinom(n, 1, 0.5) * 1.0) else matrix(0, n, 0)
dat <- lapply(seq_len(genes), function(i) {
  G <- matrix(rbinom(n * snps, 2, runif(snps, 0.05, 0.5)), n, byrow = TRUE) * 1.0
  y <- 0.3 * G[, 1] + rnorm(n)
  list(G = G, Uy = drop(crossprod(U, y)))
})
invisible(fit(U, lambda, dat[[1]]$Uy, UX0, dat[[1]]$G, C, FALSE))      # warm-up
tm <- system.time(r <- lapply(dat, function(d) fit(U, lambda, d$Uy, UX0, d$G, C, FALSE)))
cat(sprintf("suez %s  N %d  genes %d  snps %d  interaction %d  elapsed %.2f s  user %.2f s\n",
            packageVersion("suez", lib.loc = lib), n, genes, snps, inter, tm[["elapsed"]], tm[["user.self"]]))
saveRDS(list(time = tm, res = r, version = as.character(packageVersion("suez", lib.loc = lib)),
             BLAS = tryCatch(extSoftVersion()[["BLAS"]], error = function(e) NA)), out)
