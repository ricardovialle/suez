# Phase split of a 1-core run (base12): set-up (quantile_normalize, kernel_contrast_eigen) vs per-gene C++ fits.
#   Rscript bench/proto/r12/parallel/phases.R XXL XL_grm
Sys.setenv(SUEZ_HARNESS_LIB = "bench/lib/base12"); source("bench/proto/harness.R")
for (dsn in commandArgs(TRUE)) {
  d = readRDS(file.path("bench/data", paste0(dsn, ".rds")))
  K = d$sample_kernel
  if (is.null(K)) { K = outer(d$anno$individual, d$anno$individual, "=="); class(K) = "numeric" }
  t_qn = system.time(suez:::quantile_normalize(d$input))[["elapsed"]]
  t_eig = system.time(e <- suez:::kernel_contrast_eigen(K))[["elapsed"]]
  t_solver = system.time(.Call("suez_eigen_sym", K[-1, -1], PACKAGE = "suez"))[["elapsed"]]
  P = prepare(dsn)
  t_rot_x = system.time(for (g in P$genes) crossprod(P$U, g$G))[["elapsed"]]   # R BLAS U'G (no dedup), reference
  r = run_impl(baseline_fit, P)
  cat(sprintf("%-7s N=%d genes=%d snps=%d qn=%.2f kernel_eigen=%.2f (solver %.2f) fits=%.2f crossprodUG=%.2f\n",
              dsn, nrow(K), length(P$genes), r$n_snp, t_qn, t_eig, t_solver, r$time, t_rot_x))
}
