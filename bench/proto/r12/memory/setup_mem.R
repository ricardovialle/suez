# Rscript setup_mem.R <lib> <N> <zz|grm> <eigen|lapack>: phys_footprint before / peak during kernel_contrast_eigen alone
a = commandArgs(TRUE); suppressPackageStartupMessages(library(suez, lib.loc = a[1])); N = as.integer(a[2]); options(suez.eigen_solver = a[4])
fp = function() { o = system(sprintf("footprint -p %d --noCategories -f bytes", Sys.getpid()), intern = TRUE)
  as.numeric(sub(".*: ([0-9]+) B.*", "\\1", grep("phys_footprint", o, value = TRUE))) / 2^20 }
set.seed(1)
K = if (a[3] == "grm") { G = matrix(rbinom(N * 1500, 2, 0.3), N); G = scale(G); K = tcrossprod(G) / 1500; rm(G); K } else {
  id = rep(seq_len(N / 2), each = 2); Z = matrix(0, N, N / 2); Z[cbind(1:N, id)] = 1; tcrossprod(Z) }
invisible(gc()); f0 = fp()
t = system.time(e <- suez:::kernel_contrast_eigen(K))[[3]]
f1 = fp()
cat(sprintf("RES\t%s\t%d\t%s\t%s\t%.0f\t%.0f\t%.0f\t%.2f\t%.2f\n", basename(a[1]), N, a[3], a[4], f0[1], f1[2], f1[2] - f0[1], (f1[2] - f0[1]) * 2^20 / (8 * N^2), t))
