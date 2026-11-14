# Interleaved base12 vs a C++ prototype (fit_gene_proto) on the harness inputs.
# Rscript bench/proto/r12/algorithm/bench_cpp.R <proto.cpp> <dataset> <n_genes> <reps>
a = commandArgs(TRUE); ng = as.numeric(a[3]); reps = as.integer(a[4])
Sys.setenv(SUEZ_HARNESS_LIB = "bench/lib/base12"); source("bench/proto/harness.R")
Rcpp::sourceCpp(a[1], cacheDir = "bench/proto/r12/algorithm/.cache")
P = prepare(a[2], n_genes = ng)
tf = matrix(NA, reps, 2)
for (r in 1:reps) { x = run_impl(baseline_fit, P); tf[r, 1] = x$time; y = run_impl(fit_gene_proto, P); tf[r, 2] = y$time }
A = do.call(rbind, lapply(x$out, `[[`, "res")); B = do.call(rbind, lapply(y$out, `[[`, "res"))
relb = max(abs(A[, 4] - B[, 4]) / pmax(abs(A[, 4]), A[, 5]), na.rm = TRUE)
cat(sprintf("%s %s snps=%d base %.3f proto %.3f speed-up %.3fx (ratios %s) identical=%s rel_beta_tol=%.2g\n", basename(a[1]), a[2], x$n_snp,
  median(tf[, 1]), median(tf[, 2]), median(tf[, 1] / tf[, 2]), paste(round(tf[, 1] / tf[, 2], 3), collapse = " "),
  identical(A, B) && identical(sapply(x$out, `[[`, "l0"), sapply(y$out, `[[`, "l0")), relb))
print(signif(compare(x, y), 3))
