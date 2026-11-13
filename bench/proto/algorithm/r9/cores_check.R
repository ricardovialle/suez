# identical() results across 1/2/4 cores (different chunkings), the per-gene path for a non-matrix genotype
# (S3 get_cis_geno method, as the TabixFile path) vs the r8_fix build, and permutation_approach = "permute".
# Rscript cores_check.R <lib> <dataset> <n_genes>
a = commandArgs(TRUE); lib = a[1]; dsn = a[2]; ng = as.integer(a[3])
suppressPackageStartupMessages({ library(suez, lib.loc = lib); library(doMC) })
d = readRDS(file.path("bench/data", paste0(dsn, ".rds")))
inp = d$input[head(d$geneloc$geneid, ng), ]
run = function(cores, geno = d$genotype, perm = "boot", L = lib) {
  registerDoMC(cores); set.seed(42)
  r = map_interaction_qtl(inp, geno, d$geneloc, d$snploc, d$anno, sample_kernel = d$sample_kernel,
                          normalization_approach = "qq", permutation_approach = perm, cisdist = 1e5)
  as.data.frame(r) }
r1 = run(1); r2 = run(2); r4 = run(4)
cat(dsn, "rows", nrow(r1), " identical 1 vs 2 cores:", identical(r1, r2), " 1 vs 4 cores:", identical(r1, r4), "\n")
get_cis_geno.fakegeno <<- function(genotype, chrom, left, right, snploc, genome_build, ...)
  suez:::get_cis_geno.matrix(genotype$m, chrom, left, right, snploc, genome_build)
f1 = run(1, structure(list(m = d$genotype), class = "fakegeno"))
num = sapply(r1, is.numeric)
cat(dsn, "S3 (non-matrix) genotype path: same rows", identical(r1[!num], f1[!num]), " max |diff| vs chunked",
    max(abs(as.matrix(r1[num]) - as.matrix(f1[num])), na.rm = TRUE), "\n")
saveRDS(f1, file.path("bench/proto/algorithm/r9", paste0("s3path_", basename(lib), "_", dsn, ".rds")))
p1 = run(1, perm = "permute"); p4 = run(4, perm = "permute")
cat(dsn, "permute: rows", nrow(p1), " identical 1 vs 4 cores:", identical(p1, p4), "\n")
saveRDS(p1, file.path("bench/proto/algorithm/r9", paste0("permute_", basename(lib), "_", dsn, ".rds")))
