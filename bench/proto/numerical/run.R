# Rscript bench/proto/numerical/run.R <dataset> <n_genes> <reps> <proto.cpp> [more.cpp ...]
# correctness of each proto vs the installed baseline, and interleaved timing vs p0 (median ratio)
setwd("/Users/rvialle/suez"); source("bench/proto/harness.R")
a = commandArgs(TRUE); dsn = a[1]; ng = as.numeric(a[2]); reps = as.integer(a[3]); files = a[-(1:3)]
P = prepare(dsn, ng)
# "file.cpp@-O3 -march=native" compiles with extra flags (separate cache dir per spec)
load_fun = function(spec) { p = strsplit(spec, "@")[[1]]; Sys.setenv(PKG_CXXFLAGS = if (length(p) > 1) p[2] else "")
  e = new.env(); Rcpp::sourceCpp(p[1], env = e, cacheDir = file.path("bench/proto/numerical/.cache", make.names(spec)), rebuild = length(p) > 1); e$fit_gene_proto }
impl = c(list(p0 = load_fun("bench/proto/p0_baseline.cpp")), setNames(lapply(files, load_fun), basename(files)))
base = run_impl(baseline_fit, P)
for (nm in names(impl)[-1]) { cat(sprintf("== %s %s\n", dsn, nm)); print(signif(compare(base, run_impl(impl[[nm]], P)), 3)) }
tm = matrix(NA, reps, length(impl), dimnames = list(NULL, names(impl)))
for (r in seq_len(reps)) for (nm in names(impl)) tm[r, nm] = run_impl(impl[[nm]], P)$time
cat(sprintf("== %s timing (s), %d genes, %d SNPs\n", dsn, length(P$genes), base$n_snp)); print(tm)
cat("speed-up vs p0 (median of per-rep ratios):\n"); print(round(apply(tm[, 1] / tm[, -1, drop = FALSE], 2, median), 2))
