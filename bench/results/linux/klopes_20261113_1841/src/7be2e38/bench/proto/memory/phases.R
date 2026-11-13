# Memory by phase of map_interaction_qtl (step 2), replaying its set-up steps in order.
# Rscript bench/proto/memory/phases.R <lib> <dataset>
a = commandArgs(TRUE); lib = a[1]; dsn = a[2]
suppressPackageStartupMessages({ library(suez, lib.loc = lib); library(doMC) }); registerDoMC(1)
rss = function() as.numeric(system(paste("ps -o rss= -p", Sys.getpid()), intern = TRUE)) / 1024
mark = function(what) { g = gc(); cat(sprintf("%-28s RSS %6.1f MB  R heap max since last %6.1f MB  now %6.1f MB\n", what, rss(),
  sum(g[, ncol(g)]), sum(g[, 2]))); invisible(gc(reset = TRUE)) }
mark("library(suez)")
d = readRDS(file.path("bench/data", paste0(dsn, ".rds"))); mark("readRDS data")
cat(sprintf("  input %s, genotype %s, kernel %s\n", format(object.size(d$input), units = "MB"),
  format(object.size(d$genotype), units = "MB"), format(object.size(d$sample_kernel), units = "MB")))
input = suez:::quantile_normalize(d$input); mark("quantile_normalize")
anno = d$anno; K = d$sample_kernel
if (is.null(K)) { K = outer(anno$individual, anno$individual, "=="); class(K) = "numeric" }; mark("default kernel (outer)")
eig = suez:::kernel_contrast_eigen(K); mark("kernel_contrast_eigen")
set.seed(42)
res = map_interaction_qtl(d$input, d$genotype, d$geneloc, d$snploc, d$anno, sample_kernel = d$sample_kernel,
  normalization_approach = "qq", permutation_approach = "boot", cisdist = 1e5)
mark("map_interaction_qtl (whole)")
cat(sprintf("  result %d rows, %s (%.0f bytes/row)\n", nrow(res), format(object.size(res), units = "MB"), as.numeric(object.size(res)) / nrow(res)))
