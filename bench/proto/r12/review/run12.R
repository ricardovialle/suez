# One end-to-end run of an installed suez build, as bench/run.R, on a dataset of bench/data or of the reviewer's
# data dir:  Rscript bench/proto/r12/review/run12.R <lib dir> <dataset> <cores> <out.rds> [perm=boot] [n_genes=Inf] [null_seed]
#   null_seed: permute individuals' genotype labels with that seed first (null run, all associations broken).
a = commandArgs(TRUE); lib = a[1]; dsn = a[2]; cores = as.integer(a[3]); out = a[4]
perm = if (length(a) >= 5) a[5] else "boot"; ng = if (length(a) >= 6) as.numeric(a[6]) else Inf
suppressPackageStartupMessages({ library(suez, lib.loc = lib); library(doMC) })
registerDoMC(cores)
f = file.path("bench/data", paste0(dsn, ".rds")); if (!file.exists(f)) f = file.path("bench/proto/r12/review/data", paste0(dsn, ".rds"))
d = readRDS(f)
if (length(a) >= 7) { set.seed(as.integer(a[7]))
  if (!is.null(d$genotype)) colnames(d$genotype) = colnames(d$genotype)[sample(ncol(d$genotype))] }
genes = head(rownames(d$input), ng)
args = list(d$input[genes, , drop = FALSE], d$genotype, d$geneloc, d$snploc, d$anno, sample_kernel = d$sample_kernel,
            normalization_approach = "qq", permutation_approach = perm, cisdist = 1e5)
if (!is.null(d$covariates)) args$covariates = d$covariates
if (!is.null(d$local_ancestry)) args$local_ancestry = d$local_ancestry
set.seed(42)
t0 = proc.time()[["elapsed"]]
res = as.data.frame(do.call(map_interaction_qtl, args))
attr(res, "wall_s") = proc.time()[["elapsed"]] - t0
saveRDS(res, out)
cat(sprintf("RUN12\t%s\t%s\t%d\t%.3f\t%d\n", basename(lib), dsn, cores, attr(res, "wall_s"), nrow(res)))
