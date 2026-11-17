# Write a dataset of bench/data/<name>.rds in the files tensorQTL reads (run from the repository root).
#   Rscript bench/tensorqtl/export.R <name> ...      -> bench/tensorqtl/data/<name>/
# Expression is quantile-normalized per gene exactly as map_interaction_qtl(normalization_approach = "qq")
# does, so both tools see the same phenotype values.
suppressPackageStartupMessages(library(suez, lib.loc = "bench/lib/cur"))
for (nm in commandArgs(TRUE)) {
  d = readRDS(file.path("bench/data", paste0(nm, ".rds")))
  out = file.path("bench/tensorqtl/data", nm); dir.create(out, recursive = TRUE, showWarnings = FALSE)
  stopifnot(!anyNA(d$genotype))
  N = nrow(d$anno); samples = sprintf("s%05d", seq_len(N))
  qq = suez:::quantile_normalize(d$input); colnames(qq) = samples
  data.table::fwrite(data.table::data.table(gene_id = rownames(qq), qq), file.path(out, "pheno.tsv.gz"), sep = "\t")
  data.table::fwrite(data.frame(gene_id = d$geneloc$geneid, chr = d$geneloc$chr, start = d$geneloc$left,
                                end = d$geneloc$right), file.path(out, "pheno_pos.tsv"), sep = "\t")
  data.table::fwrite(data.frame(sample = samples, individual = d$anno$individual, condition = d$anno$condition),
                     file.path(out, "anno.tsv"), sep = "\t")
  if (!is.null(d$covariates))
    data.table::fwrite(data.frame(sample = samples, d$covariates), file.path(out, "covariates.tsv"), sep = "\t")
  G = d$genotype   # SNPs x individuals, int8, row-major
  con = file(file.path(out, "geno.i8"), "wb"); writeBin(as.integer(t(G)), con, size = 1); close(con)
  writeLines(colnames(G), file.path(out, "individuals.txt"))
  data.table::fwrite(data.frame(variant_id = d$snploc$snpid, chrom = d$snploc$chr, pos = d$snploc$pos),
                     file.path(out, "variants.tsv"), sep = "\t")
  data.table::fwrite(d$truth, file.path(out, "truth.tsv"), sep = "\t")
  if (!is.null(d$params$cisdist)) writeLines(format(d$params$cisdist, scientific = FALSE), file.path(out, "window.txt"))
  cat(nm, "exported\n")
}
