source("bench/proto/harness.R"); Rcpp::sourceCpp("bench/proto/round8/r8_timed.cpp")
for (ds in c("S", "M", "L", "L_grm", "M_dense", "XL")) {
  P = prepare(ds, if (ds %in% c("XL", "M_dense")) 8 else 25); trot(TRUE)
  P$U = suez:::kernel_contrast_eigen(if (is.null(readRDS(file.path("bench/data", paste0(ds, ".rds")))$sample_kernel)) { a = readRDS(file.path("bench/data", paste0(ds, ".rds")))$anno; K = outer(a$individual, a$individual, "=="); K * 1 } else readRDS(file.path("bench/data", paste0(ds, ".rds")))$sample_kernel)$vectors
  r = run_impl(fit_gene_proto, P); cat(sprintf("%-8s n=%4d total %.2fs rotation %.2fs (%.0f%%)\n", ds, length(P$lambda), r$time, trot(TRUE), 100 * trot(FALSE) / r$time))
}
