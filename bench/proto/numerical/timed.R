setwd("/Users/rvialle/suez"); source("bench/proto/harness.R")
a = commandArgs(TRUE); e = new.env(); Rcpp::sourceCpp(a[1], env = e, cacheDir = "bench/proto/numerical/.cache")
for (d in a[-1]) { P = prepare(d, as.numeric(Sys.getenv("NG", "20"))); e$timers(TRUE); t = run_impl(e$fit_gene_proto, P)$time; tt = e$timers(FALSE)
  cat(d, "total", t, " shares: rotate, data grid, data refine, boot sim, boot moments, boot grid, boot refine\n"); print(round(tt[1:7] / sum(tt), 3)) }
