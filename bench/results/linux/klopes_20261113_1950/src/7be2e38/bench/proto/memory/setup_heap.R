# R heap peak of map_interaction_qtl's N x N set-up (2 genes only): Rscript setup_heap.R <lib> <dataset>
a = commandArgs(TRUE); suppressPackageStartupMessages({ library(suez, lib.loc = a[1]); library(doMC) }); registerDoMC(1)
d = readRDS(file.path("bench/data", paste0(a[2], ".rds"))); inp = d$input[1:2, , drop = FALSE]
invisible(gc(reset = TRUE)); h0 = sum(gc()[, 2])
set.seed(42); r = map_interaction_qtl(inp, d$genotype, d$geneloc, d$snploc, d$anno, sample_kernel = d$sample_kernel)
g = gc(); fp = system(paste("footprint -p", Sys.getpid()), intern = TRUE)
cat(sprintf("%s %s N=%d: R heap before %.0f MB, max during call %.0f MB (+%.0f); %s\n", basename(a[1]), a[2], ncol(inp), h0,
  sum(g[, ncol(g)]), sum(g[, ncol(g)]) - h0, trimws(fp[grep("phys_footprint_peak", fp)])))
