# Bootstrap-distribution check over several seeds: Rscript bench/proto/r12/review/boot_seeds.R <cand lib> <dataset> <seed...>
# Per seed: KS test of the bootstrap interaction LR (2 (l_boot_interact - l_boot_geno)) and of l_boot_geno - l_geno,
# candidate vs base12; then the pooled KS p over the seeds' Fisher combination.
a = commandArgs(TRUE); cand = a[1]; dsn = a[2]; seeds = as.integer(a[-(1:2)])
d = readRDS(file.path("bench/data", paste0(dsn, ".rds")))
run = function(lib, seed) {
  r = callr_like(lib, seed); r
}
callr_like = function(lib, seed) {
  f = tempfile(fileext = ".rds")
  code = sprintf('suppressPackageStartupMessages({library(suez, lib.loc="bench/lib/%s"); library(doMC)}); registerDoMC(2);
    d = readRDS("bench/data/%s.rds"); set.seed(%d)
    args = list(d$input, d$genotype, d$geneloc, d$snploc, d$anno, sample_kernel = d$sample_kernel, normalization_approach = "qq",
                permutation_approach = "boot", cisdist = 1e5)
    if (!is.null(d$covariates)) args$covariates = d$covariates
    r = do.call(map_interaction_qtl, args); saveRDS(r[, c("l_geno", "l_boot_geno", "l_boot_interact")], "%s")', lib, dsn, seed, f)
  system2("Rscript", c("-e", shQuote(code)), stdout = FALSE, stderr = FALSE); readRDS(f)
}
ps = NULL
for (s in seeds) {
  b = callr_like("base12", s); c = callr_like(cand, s)
  lr = function(x) { v = 2 * (x$l_boot_interact - x$l_boot_geno); v[is.finite(v)] }
  dg = function(x) { v = x$l_boot_geno - x$l_geno; v[is.finite(v)] }
  k1 = suppressWarnings(ks.test(lr(b), lr(c))$p.value); k2 = suppressWarnings(ks.test(dg(b), dg(c))$p.value)
  q = function(v) round(quantile(v, c(.5, .9, .99, .999)), 3)
  cat(sprintf("%s seed %d: KS p LR %.3g, l_boot_geno-l_geno %.3g | LR q %s / %s\n", dsn, s, k1, k2,
              paste(q(lr(b)), collapse = " "), paste(q(lr(c)), collapse = " ")))
  ps = rbind(ps, c(k1, k2))
}
fisher = function(p) pchisq(-2 * sum(log(p)), 2 * length(p), lower.tail = FALSE)
cat(sprintf("%s pooled (Fisher) p: LR %.3g, l_boot_geno-l_geno %.3g\n", dsn, fisher(ps[, 1]), fisher(ps[, 2])))
