# Peak memory of one suez_fit_gene call on a wide cis window (real data: 1-7k SNPs per gene), base12.
#   /usr/bin/time -l Rscript bench/proto/r12/review/mem_gene.R <lib> <n_ind> <S>   (S = 0: set-up only)
# The difference of maximum RSS between S and S = 0 is the per-gene working set of one worker.
a = commandArgs(TRUE); suppressPackageStartupMessages(library(suez, lib.loc = a[1]))
ni = as.integer(a[2]); S = as.integer(a[3]); set.seed(1)
ind = rep(seq_len(ni), each = 2); N = length(ind); cond = rep(c(0, 1), ni)
K = outer(ind, ind, "=="); K[] = as.numeric(K)
e = suez:::kernel_contrast_eigen(K); rm(K); invisible(gc())
UX0 = crossprod(e$vectors, scale(cond)); Uy = drop(crossprod(e$vectors, rnorm(N)))
if (S > 0) {
  G = matrix(rbinom(ni * S, 2, 0.3), ni)[ind, ]; storage.mode(G) = "double"
  g0 = gc(reset = TRUE)
  t0 = proc.time()[["elapsed"]]
  f = .Call("suez_fit_gene", e$vectors, e$values, Uy, UX0, G, cbind(cond), TRUE, PACKAGE = "suez")
  cat(sprintf("N=%d S=%d fit %.2f s, G %.0f MB\n", N, S, proc.time()[["elapsed"]] - t0, 8 * N * S / 2^20))
}
