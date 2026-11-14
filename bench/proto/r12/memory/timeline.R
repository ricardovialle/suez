# Rscript timeline.R <lib> <dataset> [n_genes]: RSS (ps) and R heap peak (gc max used) after each set-up phase
a = commandArgs(TRUE); suppressPackageStartupMessages({ library(suez, lib.loc = a[1]); library(doMC) }); registerDoMC(1)
ns = asNamespace("suez")
rss = function() as.numeric(system(sprintf("ps -o rss= -p %d", Sys.getpid()), intern = TRUE)) / 1024
pt = function(s) { g = gc(); cat(sprintf("%-28s rss %7.1f MB  R heap now %7.1f  R heap max %7.1f MB\n", s, rss(), sum(g[, 2]), sum(g[, 6]))); invisible(gc(reset = TRUE)) }
pt("loaded")
d = readRDS(file.path("/Users/rvialle/suez/bench/data", paste0(a[2], ".rds"))); pt("data read")
anno = d$anno; K = d$sample_kernel
if (is.null(K)) { Z = matrix(0, nrow(anno), length(unique(anno$individual))); Z[cbind(seq_len(nrow(anno)), match(anno$individual, unique(anno$individual)))] = 1; K = tcrossprod(Z); rm(Z) }
pt("kernel built (Z Z')")
eig = ns$kernel_contrast_eigen(K); pt("kernel_contrast_eigen")
K = NULL; pt("kernel freed")
ng = if (length(a) >= 3) as.integer(a[3]) else 2; g = head(d$geneloc$geneid, ng)
set.seed(1); r = map_interaction_qtl(d$input[g, , drop = FALSE], d$genotype, d$geneloc, d$snploc, d$anno, sample_kernel = d$sample_kernel)
pt(sprintf("map_interaction_qtl %d genes", ng))
