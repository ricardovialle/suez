# Datasets for the suez vs tensorQTL comparison. Most are the suez benchmark workloads of bench/datasets.R
# (same simulator, same seeds); two are added here:
#   G     a larger scan: 400 samples (200 individuals x 2 conditions), 4,000 genes x 125 SNPs (~1M tests)
#   NULLR a null set for calibration: no cis effects, strong individual effect (ICC 0.5), 1,000 genes
source("bench/datasets.R")
TQ_DATASETS = c(DATASETS[c("S", "M", "L", "XL", "XXL", "M_grm", "L_grm", "XL_grm", "M_fac3", "M_num",
                           "M_cov", "M_eqtl", "M_unbal", "HS")],
                list(G     = ds(n_ind = 200, n_gene = 4000, snps_per_gene = 125),
                     NULLR = ds(n_ind = 200, n_gene = 1000, snps_per_gene = 50, effect = 0, ind_sd = 1),
                     # human-cohort shape: 1 gene per 200 kb, cis window +-1 Mb (each SNP in ~10 genes' windows,
                     # ~5,000 SNPs per gene, 1 per 400 bp), 500 individuals x 2 conditions, GRM kernel
                     REAL  = ds(n_ind = 500, n_gene = 150, snps_per_gene = 500, gene_spacing = 2e5, cis_span = 1e6,
                                cisdist = 1e6, kernel = "grm")))
# speed: suez default kernel unless the dataset defines one; tensorQTL treats samples as independent
SPEED = c("S", "M", "L", "XL", "XXL", "G", "M_grm", "L_grm", "XL_grm", "M_fac3", "M_num", "M_cov",
          "M_eqtl", "M_unbal")
# statistics: calibration and power against the simulated truth
STATS = c("M", "M_fac3", "M_num", "M_unbal", "M_cov", "M_eqtl", "HS", "NULLR")
