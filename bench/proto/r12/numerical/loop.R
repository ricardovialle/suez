Sys.setenv(SUEZ_HARNESS_LIB = "bench/lib/base12"); source("bench/proto/harness.R")
a = commandArgs(TRUE); env = new.env(); Rcpp::sourceCpp(a[1], env = env); P = prepare(a[2], as.integer(a[3]))
cat("READY\n"); for (i in 1:as.integer(a[4])) run_impl(env$fit_gene_proto, P)
