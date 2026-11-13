a = commandArgs(TRUE); source("bench/proto/harness.R")
Rcpp::sourceCpp(a[1], cacheDir = "bench/proto/algorithm/.cache")
for (ds in a[-(1:2)]) { P = prepare(ds, as.numeric(a[2])); timers(TRUE); r = run_impl(fit_gene_proto, P); tt = timers(FALSE)
cat(ds, sprintf("total %.3f | null+rotation %.3f  compress %.3f gridGEMM %.3f gridsolve %.3f | data(rest) %.3f sim %.3f boot(rest) %.3f\n", r$time, tt[1], tt[2], tt[3], tt[6], tt[4]-0, tt[5], tt[7])) }
