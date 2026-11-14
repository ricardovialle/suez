Sys.setenv(SUEZ_HARNESS_LIB = "bench/lib/base12"); source("bench/proto/harness.R")
Rcpp::sourceCpp("bench/proto/r12/numerical/p_instr.cpp")
for (ds in commandArgs(TRUE)) { P = prepare(ds, 10); ptimers(TRUE); r = run_impl(fit_gene_proto, P); t = ptimers(TRUE)
  cat(ds, "total", r$time, "snps", r$n_snp, "\n"); print(round(c(t[1:3], ms_per_call = 1e6 * sum(t[1:3]) / sum(t[4:5]), calls_per_refine = sum(t[4:5]) / t[6], t[4:6]), 4)) }
