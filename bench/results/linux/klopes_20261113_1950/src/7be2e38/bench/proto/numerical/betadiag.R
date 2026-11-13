# beta differences vs baseline in units of se, and relative diff where |beta| is not ~0
setwd("/Users/rvialle/suez"); source("bench/proto/harness.R")
a = commandArgs(TRUE); dsn = a[1]; P = prepare(dsn, as.numeric(a[2])); base = run_impl(baseline_fit, P)
B = do.call(rbind, lapply(base$out, `[[`, "res"))
for (f in a[-(1:2)]) { p = strsplit(f, "@")[[1]]; Sys.setenv(PKG_CXXFLAGS = if (length(p) > 1) p[2] else ""); e = new.env(); Rcpp::sourceCpp(p[1], env = e, cacheDir = file.path("bench/proto/numerical/.cache", make.names(f)), rebuild = length(p) > 1)
  A0 = run_impl(e$fit_gene_proto, P); print(signif(compare(base, A0), 3))
  A = do.call(rbind, lapply(A0$out, `[[`, "res")); d = abs(A[, 4] - B[, 4]); r = d / abs(B[, 4])
  i = which.max(r)
  cat(sprintf("%s %s: max|dbeta|/se=%.2e  max rel(|beta|>se/100)=%.2e  worst rel %.2e at beta=%.2e se=%.2e\n", dsn, basename(f),
      max(d / B[, 5], na.rm = TRUE), max(r[abs(B[, 4]) > B[, 5] / 100], na.rm = TRUE), r[i], B[i, 4], B[i, 5])) }
