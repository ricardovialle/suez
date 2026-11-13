# Simulates the example data in inst/extdata used by vignettes/suez.Rmd:
#   counts.tsv.gz   raw RNA-seq counts, genes x samples
#   samples.tsv     sample_id, individual, condition
#   genes.tsv       gene_id, chr, start, end
#   genotypes.vcf.gz (+ .tbi)  joint-called VCF, GT only
# 60 individuals x 2 conditions; 300 genes on chr22 (60 too lowly expressed to keep); ~10 SNPs per gene,
# some failing MAF or call-rate filters. 30 genes have an eQTL, 10 more a response eQTL (treated only).
# Run from the package root: Rscript data-raw/make_example_data.R
set.seed(2024)
out = "inst/extdata"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
n_ind = 60; n_gene = 300; snps_per_gene = 10
ind = sprintf("IND%03d", 1:n_ind)
samples = data.frame(sample_id = paste0(rep(ind, each = 2), "_", c("ctrl", "trt")),
                     individual = rep(ind, each = 2), condition = rep(c("control", "treated"), n_ind))
N = nrow(samples); trt = samples$condition == "treated"

genes = data.frame(gene_id = sprintf("ENSG%011d", 1:n_gene), chr = "chr22",
                   start = 2e7 + (1:n_gene) * 3e4 + sample(0:5000, n_gene))
genes$end = genes$start + sample(2000:20000, n_gene)

# SNPs around each gene; allele frequencies include rare variants
pos = sort(unique(rep(genes$start, each = snps_per_gene) + sample(-8e4:8e4, n_gene * snps_per_gene, TRUE)))
n_snp = length(pos)
af = c(runif(n_snp * 0.85, 0.05, 0.5), runif(n_snp - round(n_snp * 0.85), 0.002, 0.03))[sample(n_snp)]
G = t(sapply(af, function(p) rbinom(n_ind, 2, p)))                      # SNPs x individuals
snp_id = paste0("rs", 1e6 + seq_len(n_snp))

# expression: baseline + condition + individual + 2 hidden factors + genetic effects
lib = round(runif(N, 1.5e7, 4e7))
base = rnorm(n_gene, -11, 1.2); base[sample(n_gene, 60)] = -17          # 60 barely expressed genes
eta = base + outer(rnorm(n_gene, 0, 0.3), trt) +
  outer(rep(1, n_gene), rnorm(n_ind, 0, 0.25)[match(samples$individual, ind)]) +
  matrix(rnorm(n_gene * 2, 0, 0.3), n_gene) %*% matrix(rnorm(2 * N), 2)
nearest = function(g) { w = which(abs(pos - genes$start[g]) < 5e4 & af > 0.1); if (length(w)) w[1] else NA }
expressed = which(base > -15)
eqtl = sample(expressed, 30); reqtl = sample(setdiff(expressed, eqtl), 10)
truth = data.frame()
for (g in c(eqtl, reqtl)) {
  s = nearest(g); if (is.na(s)) next
  dose = G[s, match(samples$individual, ind)]
  eta[g, ] = eta[g, ] + if (g %in% eqtl) 0.35 * dose else 0.6 * dose * trt
  truth = rbind(truth, data.frame(gene_id = genes$gene_id[g], snp = snp_id[s],
                                  type = if (g %in% eqtl) "eQTL" else "response eQTL"))
}
mu = sweep(exp(eta), 2, lib, "*")
counts = matrix(rnbinom(length(mu), mu = mu, size = 20), n_gene,
                dimnames = list(genes$gene_id, samples$sample_id))

# VCF: phased-looking unphased GT, ~2% missing calls overall, a few SNPs with poor call rate
gt = ifelse(G == 0, "0/0", ifelse(G == 1, "0/1", "1/1"))
gt[matrix(runif(length(gt)) < 0.005, nrow(gt))] = "./."
bad = sample(n_snp, 40); gt[bad, ] = ifelse(runif(40 * n_ind) < 0.2, "./.", gt[bad, ])
ref = sample(c("A", "C", "G", "T"), n_snp, TRUE)
alt = sapply(ref, function(r) sample(setdiff(c("A", "C", "G", "T"), r), 1))
vcf_txt = file.path(tempdir(), "genotypes.vcf")
writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=chr22,length=50818468>",
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", ind), collapse = "\t"),
             paste("chr22", pos, snp_id, ref, alt, ".", "PASS", ".", "GT", apply(gt, 1, paste, collapse = "\t"), sep = "\t")),
           vcf_txt)
unlink(file.path(out, "genotypes.vcf.gz*"))
Rsamtools::indexTabix(Rsamtools::bgzip(vcf_txt, file.path(out, "genotypes.vcf.gz"), overwrite = TRUE), format = "vcf")

readr::write_tsv(cbind(data.frame(gene_id = rownames(counts)), as.data.frame(counts)), file.path(out, "counts.tsv.gz"))
readr::write_tsv(samples, file.path(out, "samples.tsv"))
readr::write_tsv(genes, file.path(out, "genes.tsv"))
readr::write_tsv(truth, file.path(out, "truth.tsv"))                     # simulated effects, for checking
cat("written:", dir(out), "\n")

# ---- admixed cohort for the local-ancestry section (inst/extdata/admixed) ----
# 150 two-way admixed individuals (AFR/EUR, ~75% AFR) x 2 conditions, 80 genes on chr22, 15 SNPs per gene,
# phased haplotypes with local-ancestry tracts and ancestry-specific allele frequencies. Written as
# Tractor's extract_tracts output (cohort.anc0/anc1 .dosage/.hapcount .txt.gz), plus counts and truth.
# Runs on its own: the files above are unchanged.
set.seed(2025)
out_a = file.path(out, "admixed"); dir.create(out_a, showWarnings = FALSE)
n_ind = 150; n_gene = 80; ns = 15
ind = sprintf("ADM%03d", 1:n_ind)
samples = data.frame(sample_id = paste0(rep(ind, each = 2), "_", c("ctrl", "trt")),
                     individual = rep(ind, each = 2), condition = rep(c("control", "treated"), n_ind))
N = nrow(samples); trt = samples$condition == "treated"; si = match(samples$individual, ind)
genes = data.frame(gene_id = sprintf("ENSG%011d", 1000 + 1:n_gene), chr = "chr22",
                   start = 2e7 + (1:n_gene) * 3e5)
genes$end = genes$start + 1e4
theta = rbeta(n_ind, 6, 2)                                             # global AFR fraction
LA = DA = DE = matrix(0L, n_gene * ns, n_ind)                          # AFR haplotypes, ALT on AFR / EUR
for (g in 1:n_gene) {
  rows = (g - 1) * ns + 1:ns
  pA = runif(ns, 0.15, 0.5); pE = pmin(pmax(pA + rnorm(ns, 0, 0.2), 0.1), 0.9)
  for (h in 1:2) {
    afr = matrix(runif(n_ind) < theta, n_ind, ns)                      # one tract per haplotype and region,
    sw = which(runif(n_ind) < 0.1)                                      # 10% with a switch inside it
    for (i in sw) { k = sample(2:ns, 1); afr[i, k:ns] = !afr[i, 1] }
    alt = matrix(runif(n_ind * ns), n_ind) < ifelse(afr, rep(pA, each = n_ind), rep(pE, each = n_ind))
    LA[rows, ] = LA[rows, ] + t(afr); DA[rows, ] = DA[rows, ] + t(alt & afr); DE[rows, ] = DE[rows, ] + t(alt & !afr)
  }
}
pos = rep(genes$start, each = ns) + rep(seq(-6e4, 6e4, length.out = ns), n_gene)
snp_id = paste0("rs", 2e6 + seq_along(pos))
lib = round(runif(N, 1.5e7, 4e7))
eta = rnorm(n_gene, -11, 1) + outer(rnorm(n_gene, 0, 0.3), trt) +
  outer(rep(1, n_gene), rnorm(n_ind, 0, 0.25)[si]) +
  outer(rnorm(n_gene, 0, 0.3), theta[si] - mean(theta)) +              # global-ancestry effect
  matrix(rnorm(n_gene * 2, 0, 0.3), n_gene) %*% matrix(rnorm(2 * N), 2)
types = rep(c("AFR eQTL", "EUR eQTL", "shared eQTL", "AFR response eQTL"), each = 8)
planted = sample(n_gene, length(types)); truth = data.frame()
for (j in seq_along(planted)) {
  g = planted[j]; s = (g - 1) * ns + 8                                  # middle SNP of the gene
  a = DA[s, si]; e = DE[s, si]
  eta[g, ] = eta[g, ] + switch(types[j], "AFR eQTL" = 0.5 * a, "EUR eQTL" = 0.7 * e,
                               "shared eQTL" = 0.4 * (a + e), "AFR response eQTL" = 0.8 * a * trt)
  truth = rbind(truth, data.frame(gene_id = genes$gene_id[g], snp = snp_id[s], type = types[j]))
}
counts = matrix(rnbinom(n_gene * N, mu = sweep(exp(eta), 2, lib, "*"), size = 20), n_gene,
                dimnames = list(genes$gene_id, samples$sample_id))
ref = sample(c("A", "C", "G", "T"), length(pos), TRUE)
alt = sapply(ref, function(r) sample(setdiff(c("A", "C", "G", "T"), r), 1))
tractor_file = function(m, f) {   # Tractor layout: CHROM POS ID REF ALT, then one column per individual
  colnames(m) = ind
  readr::write_tsv(cbind(data.frame(CHROM = "chr22", POS = pos, ID = snp_id, REF = ref, ALT = alt), as.data.frame(m)),
                   file.path(out_a, f))
}
tractor_file(DA, "cohort.anc0.dosage.txt.gz"); tractor_file(LA, "cohort.anc0.hapcount.txt.gz")
tractor_file(DE, "cohort.anc1.dosage.txt.gz"); tractor_file(2L - LA, "cohort.anc1.hapcount.txt.gz")
readr::write_tsv(cbind(data.frame(gene_id = rownames(counts)), as.data.frame(counts)), file.path(out_a, "counts.tsv.gz"))
readr::write_tsv(samples, file.path(out_a, "samples.tsv"))
readr::write_tsv(genes, file.path(out_a, "genes.tsv"))
readr::write_tsv(truth, file.path(out_a, "truth.tsv"))
cat("written:", dir(out_a), "\n")
