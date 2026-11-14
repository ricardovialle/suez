# identical() of map_interaction_qtl across options(suez.threads) 1/2/4/8 (lib = a build with suez.threads), and of
# threads = 1 against a reference build (run in its own process, cached).
#   Rscript bench/proto/r12/parallel/thr_identity.R <lib> <ref lib> <dataset:n_genes> ...
a = commandArgs(TRUE); lib = a[1]; ref = a[2]
run = function(lib, dsn, ng, nt) {
  f = tempfile(fileext = ".rds")
  code = sprintf('suppressPackageStartupMessages({library(suez, lib.loc = "bench/lib/%s"); library(doMC)}); registerDoMC(1)
    d = readRDS("bench/data/%s.rds"); g = head(intersect(rownames(d$input), d$geneloc$geneid), %d)
    args = list(d$input[g, , drop = FALSE], d$genotype, d$geneloc, d$snploc, d$anno, sample_kernel = d$sample_kernel,
                normalization_approach = "qq", permutation_approach = "boot", cisdist = 1e5)
    if (!is.null(d$covariates)) args$covariates = d$covariates
    if (!is.null(d$local_ancestry)) args$local_ancestry = d$local_ancestry
    options(suez.threads = %d); set.seed(42); t0 = proc.time()[[3]]; r = do.call(map_interaction_qtl, args)
    attr(r, "secs") = proc.time()[[3]] - t0; saveRDS(r, "%s")', lib, dsn, ng, nt, f)
  system2("Rscript", c("-e", shQuote(code)), stdout = FALSE, stderr = FALSE)
  r = readRDS(f); unlink(f); r
}
for (x in a[-(1:2)]) {
  dsn = sub(":.*", "", x); ng = as.integer(sub(".*:", "", x))
  rf = run(ref, dsn, ng, 1L); rs = lapply(c(1L, 2L, 4L, 8L), function(nt) run(lib, dsn, ng, nt))
  secs = sapply(rs, attr, "secs"); strip = function(r) { attr(r, "secs") = NULL; r }
  rf = strip(rf); rs = lapply(rs, strip)
  cat(sprintf("%-8s genes=%d rows=%d | threads 1 identical to %s: %s | 2/4/8 identical to 1: %s | secs 1/2/4/8: %s\n",
              dsn, ng, nrow(rs[[1]]), ref, identical(rs[[1]], rf),
              paste(sapply(rs[-1], identical, rs[[1]]), collapse = "/"), paste(sprintf("%.1f", secs), collapse = "/")))
}
