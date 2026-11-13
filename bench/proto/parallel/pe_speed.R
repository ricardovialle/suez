# Single-thread kernel time on 15 M genes; run once normally and once under `taskpolicy -c background`
# (clamped to efficiency cores) to get the E/P per-gene speed ratio.
suppressPackageStartupMessages({ library(suez, lib.loc = "bench/lib/baseline"); library(Rcpp) })
h = readLines("bench/proto/harness.R"); eval(parse(text = h[!grepl("^suppressPackageStartupMessages", h)]))
P = prepare("M", n_genes = 15); cat(run_impl(baseline_fit, P, reps = 3)$time, "\n")
