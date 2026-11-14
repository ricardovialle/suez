# GEMM dispatch modes of p_blas.cpp, interleaved: base12 vs mode 0 (Eigen) vs mode 2 (Accelerate) [vs mode 1 (R BLAS)]
#   Rscript cmp_modes.R <reps> <modes, e.g. 0,2> <dataset:genes> ...
Sys.setenv(SUEZ_HARNESS_LIB = "bench/lib/base12"); source("bench/proto/harness.R")
a = commandArgs(TRUE); reps = as.integer(a[1]); modes = as.integer(strsplit(a[2], ",")[[1]])
Sys.setenv(PKG_CPPFLAGS = "-DSUEZ_ACCELERATE", PKG_LIBS = "-framework Accelerate $(BLAS_LIBS)")
env = new.env(); Rcpp::sourceCpp("bench/proto/r12/numerical/p_blas.cpp", env = env, rebuild = TRUE)
for (spec in a[-(1:2)]) {
  ds = sub(":.*", "", spec); ng = as.integer(sub(".*:", "", spec)); P = prepare(ds, ng)
  T = matrix(NA, reps, 1 + length(modes)); res = list()
  for (r in seq_len(reps)) {
    b = run_impl(baseline_fit, P); T[r, 1] = b$time
    for (j in seq_along(modes)) { env$set_gemm_mode(modes[j]); p = run_impl(env$fit_gene_proto, P); T[r, j + 1] = p$time; res[[j]] = p }
  }
  cat(sprintf("%-8s N=%d snps=%d  base %.3f s; ", ds, length(P$lambda) + 1, b$n_snp, median(T[, 1])))
  for (j in seq_along(modes)) { cm = compare(b, res[[j]])
    cat(sprintf("mode %d: %.3f s ratio %.2f (dl %.1e, dl_boot %.1e, NAmis %d); ", modes[j], median(T[, j + 1]), median(T[, 1] / T[, j + 1]),
                max(cm[c("l0", "l_geno", "l_interact")]), max(cm[c("l_boot_geno", "l_boot_interact")]), cm["na_mismatch"])) }
  cat("\n")
}
