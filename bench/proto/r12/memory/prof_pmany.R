# Rscript prof_pmany.R <lib> [n_genes] [memory.profiling]: by.self profile of map_interaction_qtl on P_many
a = commandArgs(TRUE); ng = as.integer(a[2])
suppressPackageStartupMessages({ library(suez, lib.loc = a[1]); library(doMC) }); registerDoMC(1)
d = readRDS("bench/data/P_many.rds"); g = head(d$geneloc$geneid, ng)
f = tempfile(); Rprof(f, interval = 0.002, line.profiling = FALSE)
set.seed(42)
t = system.time(r <- map_interaction_qtl(d$input[g, , drop = FALSE], d$genotype, d$geneloc, d$snploc, d$anno))
Rprof(NULL); s = summaryRprof(f)
cat("wall", t[["elapsed"]], "\n"); print(head(s$by.self, 25))
