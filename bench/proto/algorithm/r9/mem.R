a = commandArgs(TRUE); suppressPackageStartupMessages({ library(suez, lib.loc = a[1]); library(doMC) }); registerDoMC(1)
d = readRDS(file.path("bench/data", paste0(a[2], ".rds"))); gc(reset = TRUE); set.seed(42)
r = map_interaction_qtl(d$input, d$genotype, d$geneloc, d$snploc, d$anno, sample_kernel = d$sample_kernel, cisdist = 1e5)
g = gc(); cat(a[1], a[2], "R heap max used MB:", sum(g[, ncol(g)]), "\n")
