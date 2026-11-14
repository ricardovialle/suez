# Peak RSS of the installed build <lib> fitting <dataset> <genes> (no compilation in the process tree)
a = commandArgs(TRUE); Sys.setenv(SUEZ_HARNESS_LIB = a[1]); source("bench/proto/harness.R")
P = prepare(a[2], as.integer(a[3])); gc(reset = TRUE); r = run_impl(baseline_fit, P); cat("time", r$time, "\n")
