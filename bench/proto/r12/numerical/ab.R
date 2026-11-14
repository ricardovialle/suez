# Interleaved A/B of two prototypes (and base12): Rscript ab.R <a.cpp> <b.cpp> <reps> <dataset:genes> ...
Sys.setenv(SUEZ_HARNESS_LIB = "bench/lib/base12"); source("bench/proto/harness.R")
a = commandArgs(TRUE); ea = new.env(); eb = new.env(); Rcpp::sourceCpp(a[1], env = ea); Rcpp::sourceCpp(a[2], env = eb); reps = as.integer(a[3])
for (spec in a[-(1:3)]) {
  ds = sub(":.*", "", spec); P = prepare(ds, as.integer(sub(".*:", "", spec))); T = matrix(0, reps, 3)
  for (r in 1:reps) { T[r, 1] = run_impl(baseline_fit, P)$time; T[r, 2] = run_impl(ea$fit_gene_proto, P)$time; T[r, 3] = run_impl(eb$fit_gene_proto, P)$time }
  cat(sprintf("%-8s base %.3f  A %.3f (x%.2f)  B %.3f (x%.2f)  B/A time %.3f\n", ds, median(T[, 1]), median(T[, 2]), median(T[, 1] / T[, 2]),
              median(T[, 3]), median(T[, 1] / T[, 3]), median(T[, 3] / T[, 2])))
}
