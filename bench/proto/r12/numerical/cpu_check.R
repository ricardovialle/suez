Sys.setenv(SUEZ_HARNESS_LIB = "bench/lib/base12"); source("bench/proto/harness.R")
Sys.setenv(PKG_CPPFLAGS = "-DSUEZ_ACCELERATE", PKG_LIBS = "-framework Accelerate $(BLAS_LIBS)")
env = new.env(); Rcpp::sourceCpp("bench/proto/r12/numerical/p_blas.cpp", env = env)
P = prepare("L_grm", 6)
for (m in c(0, 2, 2, 0)) { env$set_gemm_mode(m); t = system.time(run_impl(env$fit_gene_proto, P)); cat("mode", m, "elapsed", t[["elapsed"]], "user+sys", t[["user.self"]] + t[["sys.self"]], "\n") }
