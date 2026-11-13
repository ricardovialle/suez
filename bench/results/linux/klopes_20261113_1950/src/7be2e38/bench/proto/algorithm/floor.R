# Intrinsic sensitivity of the baseline: same code on inputs perturbed by ~1e-15 relative (rounding level).
source("bench/proto/harness.R")
for (ds in commandArgs(TRUE)) { P = prepare(ds, if (ds == "L_grm") 10 else 20); a = run_impl(baseline_fit, P)
  set.seed(3); for (i in seq_along(P$genes)) P$genes[[i]]$Uy = P$genes[[i]]$Uy * (1 + 1e-15 * rnorm(length(P$genes[[i]]$Uy)))
  b = run_impl(baseline_fit, P)
  A = do.call(rbind, lapply(a$out, `[[`, "res")); B = do.call(rbind, lapply(b$out, `[[`, "res"))
  cat(ds, paste(sprintf("%s=%.2g", names(compare(a, b)), compare(a, b)), collapse = " "), "dbeta/se=", max(abs(A[,4]-B[,4])/A[,5], na.rm=TRUE), "\n") }
