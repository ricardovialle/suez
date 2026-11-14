# Rscript edge.R <lib> <tag>: checkpoint (write + resume), permute, none, maf_threshold on S (20 genes)
a = commandArgs(TRUE); suppressPackageStartupMessages({ library(suez, lib.loc = a[1]); library(doMC) }); registerDoMC(1)
d = readRDS("/Users/rvialle/suez/bench/data/S.rds"); g = head(d$geneloc$geneid, 20); inp = d$input[g, ]
ck = file.path(tempdir(), "ck"); r = list()
run = function(...) { set.seed(1); map_interaction_qtl(inp, d$genotype, d$geneloc, d$snploc, d$anno, ...) }
r$ck1 = run(checkpoint_dir = ck); r$ck2 = run(checkpoint_dir = ck)
r$perm = run(permutation_approach = "permute"); r$none = run(permutation_approach = "none"); r$maf = run(maf_threshold = 0.1)
saveRDS(r, sprintf("out/edge_%s.rds", a[2]))
