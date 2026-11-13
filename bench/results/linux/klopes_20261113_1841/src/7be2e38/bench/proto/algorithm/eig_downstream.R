# Idea 3 downstream: same baseline fitter, U from the O(n^2) Householder path vs the dense one.
source("bench/proto/harness.R"); source("bench/proto/algorithm/eig_householder_fun.R")
for (ds in c("S", "M_grm")) {
  P = prepare(ds, 20); d = readRDS(file.path("bench/data", paste0(ds, ".rds"))); K = d$sample_kernel
  if (is.null(K)) { K = outer(d$anno$individual, d$anno$individual, "=="); class(K) = "numeric" }
  f = kce_fast(K); T = crossprod(f$vectors, P$U)
  Q = P; Q$U = f$vectors; Q$lambda = f$values; Q$UX0 = T %*% P$UX0
  for (i in seq_along(Q$genes)) Q$genes[[i]]$Uy = drop(T %*% P$genes[[i]]$Uy)
  a = run_impl(baseline_fit, P); b = run_impl(baseline_fit, Q)
  cm = compare(a, b); cat(ds, paste(sprintf("%s=%.2g", names(cm), cm), collapse = " "), "\n")
}
