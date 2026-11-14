# Accelerate GEMM in forked workers (as doMC / mclapply): identical results to serial, no hang; 2-worker scaling.
Sys.setenv(SUEZ_HARNESS_LIB = "bench/lib/base12"); source("bench/proto/harness.R")
Sys.setenv(PKG_CPPFLAGS = "-DSUEZ_ACCELERATE", PKG_LIBS = "-framework Accelerate $(BLAS_LIBS)")
env = new.env(); Rcpp::sourceCpp("bench/proto/r12/numerical/p_blas.cpp", env = env)
P = prepare("L_grm", 8)
one = function(g) { set.seed(g$seed); env$fit_gene_proto(P$U, P$lambda, g$Uy, P$UX0, g$G, P$C, P$boot) }
for (m in c(0, 2)) {
  env$set_gemm_mode(m)
  t1 = system.time(ser <- lapply(P$genes, one))[["elapsed"]]
  t = system.time(par <- parallel::mclapply(P$genes, one, mc.cores = 2))[["elapsed"]]
  cat(sprintf("mode %d: serial %.2fs, 2 workers %.2fs (speed-up %.2f), identical to serial: %s\n", m, t1, t, t1 / t, identical(ser, par)))
  ser2 = lapply(P$genes, one); cat("  run-to-run identical:", identical(ser, ser2), "\n")
}
