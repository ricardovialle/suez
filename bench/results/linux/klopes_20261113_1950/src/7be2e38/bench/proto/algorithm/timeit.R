# Rscript bench/proto/algorithm/timeit.R <cand.cpp> <n_genes> <reps> <datasets...>: compare() + interleaved timing vs p0
a = commandArgs(TRUE); source("bench/proto/harness.R"); source("bench/proto/review/interleave.R")
e0 = new.env(); Rcpp::sourceCpp("bench/proto/p0_baseline.cpp", env = e0, cacheDir = "bench/proto/algorithm/.cache")
e1 = new.env(); Rcpp::sourceCpp(a[1], env = e1, cacheDir = "bench/proto/algorithm/.cache")
for (ds in a[-(1:3)]) {
  P = prepare(ds, as.numeric(a[2])); n = length(P$lambda); mx = max(P$lambda)
  cat(sprintf("\n== %s  n=%d  distinct non-null eigenvalues=%d  null(<=1e-12 max)=%d\n", ds, n,
      length(unique(signif(P$lambda[P$lambda > 1e-12 * mx], 11))), sum(P$lambda <= 1e-12 * mx)))
  b = run_impl(baseline_fit, P); e1$proto_stats(TRUE); c = run_impl(e1$fit_gene_proto, P)
  A = do.call(rbind, lapply(b$out, `[[`, "res")); B = do.call(rbind, lapply(c$out, `[[`, "res"))
  cm = compare(b, c); cm = c(cm, dbeta_over_maxbeta_se = max(abs(A[,4]-B[,4]) / pmax(abs(A[,4]), A[,5]), na.rm = TRUE))
  cat(paste(sprintf("%s=%.2g", names(cm), cm), collapse = " "), " | evals/exact/gridfallback/griddirect", paste(e1$proto_stats(FALSE), collapse = "/"), "\n")
  time_interleaved(list(p0 = e0$fit_gene_proto, cand = e1$fit_gene_proto), P, reps = as.integer(a[3]))
}
