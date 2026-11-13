# Rscript tint.R <dataset> <n_genes> <reps> cand.cpp...   correctness (compare vs installed baseline) + interleave.R timing vs sourceCpp p0
setwd("/Users/rvialle/suez"); source("bench/proto/harness.R"); source("bench/proto/review/interleave.R")
a = commandArgs(TRUE); P = prepare(a[1], as.numeric(a[2]))
l = sort(P$lambda); g = sum(c(TRUE, diff(l) > 1e-12 * max(abs(l))))
cat(sprintf("== %s: n=%d, eigenvalue groups=%d, SNPs=%d\n", a[1], length(l), g, sum(sapply(P$genes, function(x) ncol(x$G)))))
ld = function(f) { e = new.env(); Rcpp::sourceCpp(f, env = e, cacheDir = "bench/proto/numerical/.cache"); e$fit_gene_proto }
impl = c(list(p0 = ld("bench/proto/p0_baseline.cpp")), setNames(lapply(a[-(1:3)], ld), basename(a[-(1:3)])))
base = run_impl(baseline_fit, P)
for (nm in names(impl)[-1]) { cat(nm, ":\n"); print(signif(compare(base, run_impl(impl[[nm]], P)), 3)) }
time_interleaved(impl, P, reps = as.integer(a[3]))
