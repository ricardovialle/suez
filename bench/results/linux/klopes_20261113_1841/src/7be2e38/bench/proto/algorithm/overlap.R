# Idea 5: how often is the same SNP rotated (in how many genes' cis windows) per dataset?
for (dsn in c("S", "M", "L", "XL", "M_dense", "M_grm")) {
  d = readRDS(file.path("bench/data", paste0(dsn, ".rds")))
  cnt = 0; used = integer(0)
  for (i in seq_len(nrow(d$geneloc))) {
    gl = d$geneloc[i, ]
    sel = which(d$snploc$chr == gl$chr & d$snploc$pos > gl$left - 1e5 & d$snploc$pos < gl$right + 1e5)
    cnt = cnt + length(sel); used = union(used, sel)
  }
  cat(sprintf("%-8s N=%d  SNP-gene pairs %d  distinct SNPs %d  overlap factor %.2f\n", dsn, nrow(d$anno), cnt, length(used), cnt / length(used)))
}
