# Datasets for the suez vs tensorQTL comparison. Most are the suez benchmark workloads of bench/datasets.R
# (same simulator, same seeds); two are added here:
#   G     a larger scan: 400 samples (200 individuals x 2 conditions), 4,000 genes x 125 SNPs (~1M tests)
#   NULLR a null set for calibration: no cis effects, strong individual effect (ICC 0.5), 1,000 genes
source("bench/datasets.R")
TQ_DATASETS = c(DATASETS[c("S", "M", "L", "XL", "XXL", "M_grm", "L_grm", "XL_grm", "M_fac3", "M_num",
                           "M_cov", "M_eqtl", "M_unbal", "HS")],
                list(G     = ds(n_ind = 200, n_gene = 4000, snps_per_gene = 125),
                     NULLR = ds(n_ind = 200, n_gene = 1000, snps_per_gene = 50, effect = 0, ind_sd = 1)))
# speed: suez default kernel unless the dataset defines one; tensorQTL treats samples as independent
SPEED = c("S", "M", "L", "XL", "XXL", "G", "M_grm", "L_grm", "XL_grm", "M_fac3", "M_num", "M_cov",
          "M_eqtl", "M_unbal")
# statistics: calibration and power against the simulated truth
STATS = c("M", "M_fac3", "M_num", "M_unbal", "M_cov", "M_eqtl", "HS", "NULLR")
