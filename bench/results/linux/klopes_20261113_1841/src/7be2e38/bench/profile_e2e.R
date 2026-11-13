# R-level profile of map_interaction_qtl for a build: Rscript bench/profile_e2e.R <lib> <dataset> [n_genes]
a = commandArgs(TRUE); ng = if (length(a) >= 3) as.integer(a[3]) else Inf
suppressPackageStartupMessages({ library(suez, lib.loc = file.path("bench/lib", a[1])); library(doMC) }); registerDoMC(1)
d = readRDS(file.path("bench/data", paste0(a[2], ".rds")))
g = head(d$geneloc$geneid, ng)
f = tempfile(); Rprof(f, interval = 0.005)
t = system.time(do.call(map_interaction_qtl, c(list(d$input[g, , drop = FALSE], d$genotype, d$geneloc, d$snploc, d$anno,
  sample_kernel = d$sample_kernel), if (!is.null(d$covariates)) list(covariates = d$covariates))))
Rprof(NULL); s = summaryRprof(f)
cat(a[2], "wall", t[["elapsed"]], "\n"); print(head(s$by.total[, 1:2], 25))
