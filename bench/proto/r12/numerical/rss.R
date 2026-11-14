# Peak RSS of one process fitting a dataset with base12 or the prototype (mode 0 Eigen / 2 Accelerate):
#   /usr/bin/time -l Rscript rss.R <base|0|2> <dataset> <genes>
Sys.setenv(SUEZ_HARNESS_LIB = "bench/lib/base12"); source("bench/proto/harness.R")
a = commandArgs(TRUE); P = prepare(a[2], as.integer(a[3]))
if (a[1] == "base") f = baseline_fit else {
  Sys.setenv(PKG_CPPFLAGS = "-DSUEZ_ACCELERATE", PKG_LIBS = "-framework Accelerate $(BLAS_LIBS)")
  env = new.env(); Rcpp::sourceCpp("bench/proto/r12/numerical/p_blas.cpp", env = env); env$set_gemm_mode(as.integer(a[1])); f = env$fit_gene_proto }
gc(reset = TRUE); r = run_impl(f, P); cat("time", r$time, "\n")
