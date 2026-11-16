# Every code example on doc/usage.html, run end to end by `Rscript doc/src/examples.R` (about a minute).
# doc/src/build.py copies each `## @ex <name>` ... `## @end` block into the page; `## @hidden` blocks run but are
# not shown (they make the simulated studies the examples use). Keep the examples runnable: edit them here.

## @hidden setup
suppressPackageStartupMessages({ library(suez); library(doMC); library(edgeR); library(VariantAnnotation) })
t_start = Sys.time()
show = function(x) print(utils::head(x, 3))
## @end

## @ex install
# install.packages("pak")
# pak::pak("ricardovialle/suez")
# pak::pak(c("VariantAnnotation", "GenomicRanges", "IRanges"))   # only to read VCFs
library(suez)
library(doMC); registerDoMC(4)            # genes run in parallel on 4 forked workers
## @end

## @ex inputs
ext = function(f) system.file("extdata", f, package = "suez")   # bundled example study
counts  = read.delim(ext("counts.tsv.gz"), row.names = 1, check.names = FALSE)   # genes x samples
samples = read.delim(ext("samples.tsv"))       # sample_id, individual, condition
genes   = read.delim(ext("genes.tsv"))         # gene_id, chr, start, end
stopifnot(identical(colnames(counts), samples$sample_id))   # suez matches samples by position

library(edgeR)                                 # TMM-normalized log2 CPM of expressed genes
dge = DGEList(counts, group = samples$condition)
dge = dge[filterByExpr(dge), , keep.lib.sizes = FALSE]
logcpm = cpm(normLibSizes(dge), log = TRUE)
## @end

## @ex genotypes
library(VariantAnnotation)
vcf = readVcf(ext("genotypes.vcf.gz"))
vcf = vcf[isSNV(vcf)]                          # biallelic SNPs
gt = geno(vcf)$GT
dosage = matrix(c("0/0" = 0, "0|0" = 0, "0/1" = 1, "0|1" = 1, "1|0" = 1, "1/1" = 2, "1|1" = 2)[gt],
                nrow(gt), dimnames = dimnames(gt))   # SNPs x individuals, ALT allele counts, NA kept
af = rowMeans(dosage, na.rm = TRUE) / 2
dosage = dosage[rowMeans(!is.na(dosage)) >= 0.95 & pmin(af, 1 - af) >= 0.05, ]

snploc  = data.frame(snpid = rownames(dosage),
                     chr = as.character(seqnames(rowRanges(vcf)[rownames(dosage)])),
                     pos = start(rowRanges(vcf)[rownames(dosage)]))
geneloc = data.frame(geneid = genes$gene_id, chr = genes$chr, left = genes$start, right = genes$end)
## @end

## @ex step1
anno = data.frame(individual = samples$individual, condition = samples$condition)
qq = t(apply(logcpm, 1, function(x) qnorm((rank(x) - 0.5) / length(x))))   # normal scores per gene

Z = scale(t(dosage)); Z[is.na(Z)] = 0          # GRM; in a real study: genome-wide, LD-pruned SNPs
grm = tcrossprod(Z) / ncol(Z)

set.seed(1)
step1 = learn_covariance(qq, anno, kinship_matrix = grm, P = 10)
step1$optimized_fit$return_code                # 0: the optimizer converged
round(step1$optimized_fit$par$s, 3)            # noise, same-individual and kinship variances
## @end
stopifnot(step1$optimized_fit$return_code == 0)

## @ex two_conditions
set.seed(1)
res = map_interaction_qtl(logcpm, dosage, geneloc, snploc, anno,
                          sample_kernel = step1$cov_suez, cisdist = 1e5)
head(res[order(res$p_interact), c("gene", "cis_snp", "p_geno", "p_interact", "beta_geno", "se_geno")])
## @end
stopifnot(nrow(res) > 1000, min(res$p_interact) < 1e-6)

## @ex vcf_direct
res_vcf = map_interaction_qtl(logcpm, Rsamtools::TabixFile(ext("genotypes.vcf.gz")), geneloc, snploc, anno,
                              sample_kernel = step1$cov_suez, maf_threshold = 0.05)
## @end
stopifnot(nrow(res_vcf) > 1000)

## @ex eqtl_only
ctrl = samples$condition == "control"          # one sample per individual
anno_ctrl = data.frame(individual = samples$individual[ctrl], condition = 1)   # one value: eQTL tests only
res_ctrl = map_interaction_qtl(logcpm[, ctrl], dosage, geneloc, snploc, anno_ctrl,
                               sample_kernel = grm[anno_ctrl$individual, anno_ctrl$individual])
summary(res_ctrl$p_interact)                   # all NA: no interaction terms
## @end
stopifnot(all(is.na(res_ctrl$p_interact)), nrow(res_ctrl) > 1000)

## @ex multiple_testing
gene_level = function(res, p) {                # Bonferroni within gene, then BH across genes
  top = do.call(rbind, lapply(split(res, res$gene), function(d) {
    i = which.min(d[[p]])
    data.frame(gene = d$gene[i], snp = d$cis_snp[i], p = d[[p]][i], p_gene = min(1, d[[p]][i] * nrow(d)))
  }))
  top$fdr = p.adjust(top$p_gene, "BH")
  top[order(top$p), ]
}
egenes  = gene_level(res, "p_geno")
regenes = gene_level(res, "p_interact")
c(eGenes = sum(egenes$fdr < 0.05), response_eGenes = sum(regenes$fdr < 0.05))
## @end

## @ex calibration
lambda = function(p) median(qchisq(p, 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, 1)
c(geno = lambda(res$p_geno), interact = lambda(res$p_interact))   # close to 1 when calibrated

lr      = 2 * (res$l_interact - res$l_geno)              # observed interaction likelihood ratios
lr_null = 2 * (res$l_boot_interact - res$l_boot_geno)    # one bootstrap null draw per SNP
res$p_interact_boot = sapply(lr, function(x) mean(c(lr_null, Inf) >= x, na.rm = TRUE))
## @end

## @ex permute
set.seed(2)
res_null = map_interaction_qtl(logcpm, dosage, geneloc, snploc, anno, sample_kernel = step1$cov_suez,
                               permutation_approach = "permute")   # genotypes shuffled per gene
c(geno = lambda(res_null$p_geno), interact = lambda(res_null$p_interact))
## @end

## @hidden toy study for the designs the bundled data does not have
set.seed(3)
n_ind = 120; n_gene = 40; n_snp = 8
ids = sprintf("D%03d", 1:n_ind)
meta = data.frame(individual = rep(ids, each = 3), timepoint = rep(c("0h", "6h", "24h"), n_ind))
meta$dose = c("0h" = 0, "6h" = 1, "24h" = 10)[meta$timepoint] * runif(nrow(meta), 0.5, 1.5)
age0 = runif(n_ind, 20, 85); meta$age = age0[match(meta$individual, ids)] + c("0h" = 0, "6h" = 3, "24h" = 7)[meta$timepoint]
meta$sex = sample(c("F", "M"), n_ind, TRUE)[match(meta$individual, ids)]
meta$batch = sample(c("b1", "b2", "b3"), nrow(meta), TRUE)
geno = matrix(rbinom(n_gene * n_snp * n_ind, 2, 0.3), n_gene * n_snp, n_ind,
              dimnames = list(paste0("s", 1:(n_gene * n_snp)), ids))
snploc_t  = data.frame(snpid = rownames(geno), chr = "chr2", pos = rep(1:n_gene, each = n_snp) * 1e6 + 1:n_snp * 100)
geneloc_t = data.frame(geneid = paste0("G", 1:n_gene), chr = "chr2", left = 1:n_gene * 1e6, right = 1:n_gene * 1e6 + 5e3)
u = rnorm(n_ind, 0, 0.6)[match(meta$individual, ids)]
expr = t(sapply(1:n_gene, function(g) {
  gs = geno[(g - 1) * n_snp + 1, meta$individual]
  eff = switch(g %% 4 + 1, 0.4 * gs, 0.5 * gs * (meta$timepoint != "0h"), 0.3 * gs * log10(meta$dose + 1),
               0.7 * gs * exp(-((meta$age - 55) / 12)^2))
  8 + 0.5 * (meta$timepoint == "24h") + c(b1 = 0, b2 = 0.8, b3 = -0.8)[meta$batch] + eff + u + rnorm(nrow(meta), 0, 0.8)
}))
dimnames(expr) = list(geneloc_t$geneid, paste(meta$individual, meta$timepoint, sep = "_"))
## @end

## @ex factor
# meta: one row per sample (individual, timepoint, dose, age, sex, batch), in the column order of expr
meta$timepoint = factor(meta$timepoint, levels = c("0h", "6h", "24h"))   # first level is the reference
anno_tp = data.frame(individual = meta$individual, condition = meta$timepoint)
res_tp = map_interaction_qtl(expr, geno, geneloc_t, snploc_t, anno_tp)
table(res_tp$df_interact)                      # 2 df: 6h and 24h each against 0h
## @end
stopifnot(all(res_tp$df_interact == 2))

## @ex numeric
anno_dose = data.frame(individual = meta$individual, condition = log10(meta$dose + 1))
res_dose = map_interaction_qtl(expr, geno, geneloc_t, snploc_t, anno_dose)   # 1 df: linear in log dose
## @end
stopifnot(all(res_dose$df_interact == 1))

## @ex spline
B = splines::ns(meta$age, df = 3)              # natural-spline basis of age, 3 columns
anno_age = data.frame(individual = meta$individual)
anno_age$condition = B
ages = seq(30, 85, 5)
at = predict(B, ages); rownames(at) = ages     # evaluate the same basis on a grid of ages
res_age = map_interaction_qtl(expr, geno, geneloc_t, snploc_t, anno_age, curve_at = at)
hit = res_age[which.min(res_age$p_interact), ]
round(rbind(beta = unlist(hit[paste0("beta_curve_", ages)]),
            se   = unlist(hit[paste0("se_beta_curve_", ages)])), 2)
## @end
stopifnot(all(res_age$df_interact == 3), "mu_curve_50" %in% names(res_age))

## @ex unbalanced
keep = sort(sample(nrow(meta), 0.8 * nrow(meta)))   # e.g. failed libraries: 1 to 3 samples per individual
table(table(meta$individual[keep]))
anno_u = data.frame(individual = meta$individual[keep], condition = meta$timepoint[keep])
res_u = map_interaction_qtl(expr[, keep], geno, geneloc_t, snploc_t, anno_u)
## @end
stopifnot(nrow(res_u) > 0)

## @ex covariates
covs = data.frame(sex = meta$sex, batch = meta$batch)   # factors are expanded with model.matrix
res_cov = map_interaction_qtl(expr, geno, geneloc_t, snploc_t, anno_tp, covariates = covs)
## @end
stopifnot(nrow(res_cov) > 0)

## @ex step1_toy
qq_t = t(apply(expr, 1, function(x) qnorm((rank(x) - 0.5) / length(x))))
Zt = scale(t(geno)); grm_t = tcrossprod(Zt) / ncol(Zt)
set.seed(1)
cov_tp = learn_covariance(qq_t, anno_tp, kinship_matrix = grm_t, P = 5)
res_tp_k = map_interaction_qtl(expr, geno, geneloc_t, snploc_t, anno_tp, sample_kernel = cov_tp$cov_suez,
                               covariates = covs)
## @end
stopifnot(nrow(res_tp_k) > 0)

## @ex grm_kernel
K = grm_t[meta$individual, meta$individual]    # individuals x individuals GRM, expanded to samples
res_grm = map_interaction_qtl(expr, geno, geneloc_t, snploc_t, anno_tp, sample_kernel = K)
## @end
stopifnot(nrow(res_grm) > 0)

## @ex local_ancestry
adm = function(f) ext(file.path("admixed", f))
la = read_tractor(sub("\\.anc0\\.dosage\\.txt\\.gz$", "", adm("cohort.anc0.dosage.txt.gz")),
                  ancestry_names = c("AFR", "EUR"))   # dosage, hapcount, snploc
samples_a = read.delim(adm("samples.tsv"))
genes_a   = read.delim(adm("genes.tsv"))
counts_a  = read.delim(adm("counts.tsv.gz"), row.names = 1, check.names = FALSE)
dge_a = DGEList(counts_a, group = samples_a$condition)
logcpm_a = cpm(normLibSizes(dge_a[filterByExpr(dge_a), , keep.lib.sizes = FALSE]), log = TRUE)

Za = scale(t(la$dosage$AFR + la$dosage$EUR)); grm_a = tcrossprod(Za) / ncol(Za)
global = colMeans(la$hapcount$AFR) / 2         # AFR fraction per individual
anno_a = data.frame(individual = samples_a$individual, condition = samples_a$condition)
qq_a = t(apply(logcpm_a, 1, function(x) qnorm((rank(x) - 0.5) / length(x))))
set.seed(1)
step1_a = learn_covariance(qq_a, anno_a, kinship_matrix = grm_a, P = 5)

geneloc_a = data.frame(geneid = genes_a$gene_id, chr = genes_a$chr, left = genes_a$start, right = genes_a$end)
res_a = map_interaction_qtl(logcpm_a, NULL, geneloc_a, la$snploc, anno_a,
                            sample_kernel = step1_a$cov_suez, local_ancestry = la,
                            covariates = cbind(AFR = global[anno_a$individual]))
head(res_a[order(res_a$p_geno), c("gene", "cis_snp", "p_geno", "p_het", "beta_AFR", "beta_EUR")])
## @end
stopifnot(all(c("p_het", "beta_AFR", "beta_interact_EUR") %in% names(res_a)))

## @ex local_ancestry_pooled
res_pooled = map_interaction_qtl(logcpm_a, la$dosage$AFR + la$dosage$EUR, geneloc_a, la$snploc, anno_a,
                                 sample_kernel = step1_a$cov_suez,
                                 covariates = cbind(AFR = global[anno_a$individual]))
## @end
stopifnot(nrow(res_pooled) > 0)

## @ex scaling
library(doMC); registerDoMC(8)                 # 8 forked workers (doParallel on Windows)
options(suez.threads = 2)                      # and 2 threads per gene: 16 CPUs in all
set.seed(1)                                    # same results on any number of workers or threads
res_par = map_interaction_qtl(logcpm, dosage, geneloc, snploc, anno, sample_kernel = step1$cov_suez,
                              permutation_approach = "none",            # skip the bootstrap: ~1/3 faster
                              checkpoint_dir = file.path(tempdir(), "suez_ckpt"))   # resumable
options(suez.threads = 1)
## @end

## @ex by_chromosome
for (chr in unique(geneloc$chr)) {             # or one cluster job per chromosome
  out = map_interaction_qtl(logcpm, dosage, geneloc[geneloc$chr == chr, ], snploc, anno,
                            sample_kernel = step1$cov_suez)
  readr::write_tsv(out, file.path(tempdir(), paste0("suez_", chr, ".tsv.gz")))
}
## @end

cat("all examples ran in", round(as.numeric(difftime(Sys.time(), t_start, units = "secs"))), "s\n")
