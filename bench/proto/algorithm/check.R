# Rscript bench/proto/algorithm/check.R <proto.cpp> <n_genes> <reps> <datasets...>
a = commandArgs(TRUE); src = a[1]; ng = as.numeric(a[2]); reps = as.integer(a[3]); dss = a[-(1:3)]
source("bench/proto/harness.R")
Rcpp::sourceCpp(src, rebuild = FALSE, cacheDir = "bench/proto/algorithm/.cache")
proto = fit_gene_proto
p0env = new.env(); Rcpp::sourceCpp("bench/proto/p0_baseline.cpp", env = p0env, cacheDir = "bench/proto/algorithm/.cache")
for (ds in dss) {
  P = prepare(ds, ng)
  tb = tp = numeric(0)
  for (r in seq_len(reps)) {   # interleaved
    b = run_impl(p0env$fit_gene_proto, P); tb = c(tb, b$time)
    if (exists("proto_stats")) proto_stats(TRUE)
    p = run_impl(proto, P); tp = c(tp, p$time)
  }
  cmp = compare(b, p)
  A = do.call(rbind, lapply(b$out, `[[`, "res")); B = do.call(rbind, lapply(p$out, `[[`, "res"))
  cmp = c(cmp, dbeta_over_se = max(abs(A[, 4] - B[, 4]) / A[, 5], na.rm = TRUE))
  st = if (exists("proto_stats")) proto_stats(FALSE) else NA
  cat(sprintf("%-8s genes=%d snps=%d  p0 %.3fs  proto %.3fs  ratio %.2fx | %s | stats %s\n", ds, length(P$genes), b$n_snp,
              median(tb), median(tp), median(tb) / median(tp),
              paste(sprintf("%s=%.2g", names(cmp), cmp), collapse = " "), paste(st, collapse = "/")))
}
