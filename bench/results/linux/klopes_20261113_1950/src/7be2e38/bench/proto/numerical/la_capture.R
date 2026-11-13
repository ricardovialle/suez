# capture the per-gene inputs of suez_fit_models from map_interaction_qtl(local_ancestry = ...) on LA (first genes)
suppressPackageStartupMessages(library(suez, lib.loc = "bench/lib/baseline"))
d = readRDS("bench/data/LA.rds"); ng = 15
cap = new.env(); cap$x = list()
trace("la_gene", where = asNamespace("suez"), print = FALSE, exit = quote({
  cap$x[[length(cap$x) + 1]] <- list(U = U, lambda = lambda, Uy = Uy, UX0 = UX0, cols = cols, models = models, boot = do_boot) }))
geneloc = head(d$geneloc, ng)
set.seed(42)
res = map_interaction_qtl(d$input[geneloc$geneid, ], d$genotype, geneloc, d$snploc, d$anno, sample_kernel = d$sample_kernel,
        normalization_approach = "qq", permutation_approach = "boot", cisdist = 1e5, local_ancestry = d$local_ancestry)
cat("captured", length(cap$x), "genes; SNPs:", sapply(cap$x, function(z) ncol(z$cols[[1]])), "\n")
saveRDS(cap$x, "bench/proto/numerical/la_inputs.rds")
