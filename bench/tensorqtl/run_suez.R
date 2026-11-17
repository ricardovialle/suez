# Time suez on one dataset and save per-pair results (run from the repository root).
#   Rscript bench/tensorqtl/run_suez.R <dataset> <cores> <kernel> <perm> <out_prefix>
#     kernel: "dataset" (the dataset's sample_kernel; NULL = repeated-individual kernel only)
#             "learned" (step 1: learn_covariance(qq, anno, P = 10), timed separately)
#     perm:   "none" or "boot" (permutation_approach)
# cis window: the dataset's params$cisdist (default 1e5); env SUEZ_ROTATE_MEM sets options(suez.rotate_mem)
# Writes <out_prefix>.tsv.gz (gene, snp, p_geno, p_interact, beta_geno, se_geno) and <out_prefix>.json (timings).
a = commandArgs(TRUE); dsn = a[1]; cores = as.integer(a[2]); kernel = a[3]; perm = a[4]; out = a[5]
suppressPackageStartupMessages({ library(suez, lib.loc = "bench/lib/cur"); library(doMC) })
# monotonic clock: proc.time()[["elapsed"]] follows the wall clock, which jumps under WSL2
now = function() as.numeric(strsplit(readLines("/proc/uptime"), " ")[[1]][1])
registerDoMC(cores)
m = Sys.getenv("SUEZ_ROTATE_MEM"); if (nzchar(m)) options(suez.rotate_mem = as.numeric(m))
d = readRDS(file.path("bench/data", paste0(dsn, ".rds")))
cisdist = if (is.null(d$params$cisdist)) 1e5 else d$params$cisdist
t_step1 = NA
K = d$sample_kernel
if (kernel == "learned") {
  set.seed(7)
  t0 = now()
  s1 = learn_covariance(suez:::quantile_normalize(d$input), d$anno, P = 10)
  t_step1 = now() - t0
  K = s1$cov_suez
}
gc(reset = TRUE)
set.seed(42)
args = list(d$input, d$genotype, d$geneloc, d$snploc, d$anno, sample_kernel = K,
            normalization_approach = "qq", permutation_approach = perm, cisdist = cisdist)
if (!is.null(d$covariates)) args$covariates = d$covariates
t0 = now()
res = do.call(map_interaction_qtl, args)
t_step2 = now() - t0
keep = intersect(c("gene", "cis_snp", "p_geno", "p_interact", "beta_geno", "se_geno"), names(res))
data.table::fwrite(res[, keep], paste0(out, ".tsv.gz"), sep = "\t")
cat(jsonlite::toJSON(list(tool = "suez", dataset = dsn, cores = cores, kernel = kernel, perm = perm,
                          t_step1 = t_step1, t_step2 = t_step2, n_tests = nrow(res),
                          n_samples = nrow(d$anno)), auto_unbox = TRUE, digits = 6),
    file = paste0(out, ".json"))
cat(sprintf("suez %s cores=%d kernel=%s perm=%s: step1 %.2f s, step2 %.2f s, %d tests\n",
            dsn, cores, kernel, perm, t_step1, t_step2, nrow(res)))
