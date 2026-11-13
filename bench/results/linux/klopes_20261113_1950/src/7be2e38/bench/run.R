# One benchmark run: Rscript bench/run.R <lib> <dataset> <cores> <out.rds> [step]
#   lib: library holding the suez build to test (e.g. bench/lib/baseline); dataset: name in bench/data;
#   step: "step2" (map_interaction_qtl, default) or "step1" (learn_covariance).
# Prints one tab-separated line: wall, user, sys (s, this process and children) and the number of tests.
# Peak RAM and I/O come from the wrapper (bench/bench.sh runs this under /usr/bin/time -l).
a = commandArgs(TRUE); lib = a[1]; dsn = a[2]; cores = as.integer(a[3]); out = a[4]
step = if (length(a) >= 5) a[5] else "step2"
suppressPackageStartupMessages({ library(suez, lib.loc = lib); library(doMC) })
registerDoMC(cores)
d = readRDS(file.path("bench/data", paste0(dsn, ".rds")))
gc(reset = TRUE)
set.seed(42)
t0 = proc.time()
if (step == "step2") {
  args = list(d$input, d$genotype, d$geneloc, d$snploc, d$anno, sample_kernel = d$sample_kernel,
              normalization_approach = "qq", permutation_approach = "boot", cisdist = 1e5)
  if (!is.null(d$covariates)) args$covariates = d$covariates
  if (!is.null(d$local_ancestry)) args$local_ancestry = d$local_ancestry
  res = do.call(map_interaction_qtl, args)
  n = nrow(res)
} else {
  set.seed(42)
  res = learn_covariance(suez:::quantile_normalize(d$input), d$anno, P = 10)
  n = 1
}
t1 = proc.time() - t0
saveRDS(res, out)
cat(sprintf("RESULT\t%.3f\t%.3f\t%.3f\t%d\n", t1[["elapsed"]], t1[["user.self"]] + t1[["user.child"]],
            t1[["sys.self"]] + t1[["sys.child"]], n))
