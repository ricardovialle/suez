# Benchmark workloads. Defaults below; each entry overrides some of them.
base = list(seed = 1, n_ind = 100, reps = 2, condition = "factor2", n_gene = 200, snps_per_gene = 50,
            gene_spacing = 1e5, cis_span = 9e4, af_min = 0.05, ld_rho = 0.8, miss = 0, n_factor = 5,
            n_cov = 0, kernel = "default", effect = 0.4, unbalanced = FALSE)
ds = function(...) modifyList(base, list(...))
DATASETS = list(
  # size ladder (default kernel, 2-level factor, bootstrap)
  S        = ds(),                                                    # N = 200, 10k tests
  M        = ds(n_ind = 200, n_gene = 300, snps_per_gene = 100),      # N = 400, 30k tests
  L        = ds(n_ind = 500, n_gene = 150, snps_per_gene = 200),      # N = 1000, 30k tests
  XL       = ds(n_ind = 1000, n_gene = 60, snps_per_gene = 200),      # N = 2000, 12k tests
  # characteristics (N = 400 unless noted)
  M_grm    = ds(n_ind = 200, n_gene = 150, snps_per_gene = 100, kernel = "grm"),
  M_learn  = ds(n_ind = 200, n_gene = 150, snps_per_gene = 100, kernel = "learned"),
  M_num    = ds(n_ind = 200, n_gene = 150, snps_per_gene = 100, condition = "numeric", kernel = "grm"),
  M_fac3   = ds(n_ind = 200, n_gene = 150, snps_per_gene = 100, condition = "factor3", reps = 3),
  M_eqtl   = ds(n_ind = 200, n_gene = 150, snps_per_gene = 100, condition = "none", reps = 1, kernel = "grm"),
  M_cov    = ds(n_ind = 200, n_gene = 150, snps_per_gene = 100, n_cov = 10),
  M_miss   = ds(n_ind = 200, n_gene = 150, snps_per_gene = 100, miss = 0.02),
  M_dense  = ds(n_ind = 200, n_gene = 150, snps_per_gene = 100, gene_spacing = 2e4),  # overlapping cis windows
  S_rare   = ds(af_min = 0.005, ld_rho = 0.95),
  # added after the reviewer's audit (H4, H6): non-default kernels at large N, unbalanced repeats,
  # a kernel learned by learn_covariance, and a high-signal set for statistical validation
  L_grm    = ds(n_ind = 500, n_gene = 100, snps_per_gene = 200, kernel = "grm"),
  M_unbal  = ds(n_ind = 200, n_gene = 150, snps_per_gene = 100, unbalanced = TRUE),
  M_lc     = ds(n_ind = 200, n_gene = 150, snps_per_gene = 100, kernel = "lc"),
  HS       = ds(n_ind = 200, n_gene = 200, snps_per_gene = 50, effect = 1.0, kernel = "grm"),
  # local ancestry (2 ancestries)
  LA       = ds(type = "la", n_ind = 150, n_gene = 100, snps_per_gene = 50),
  # round 12: large N with a GRM kernel, very large N, many covariates, many small genes
  XL_grm   = ds(n_ind = 1000, n_gene = 40, snps_per_gene = 200, kernel = "grm"),        # N = 2000
  XXL      = ds(n_ind = 2000, n_gene = 16, snps_per_gene = 200),                        # N = 4000
  M_cov30  = ds(n_ind = 200, n_gene = 100, snps_per_gene = 100, n_cov = 30),
  P_many   = ds(n_ind = 100, n_gene = 3000, snps_per_gene = 10, gene_spacing = 3e4, cis_span = 1e4)
)
