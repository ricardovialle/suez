# Round-12 reviewer datasets (gaps of bench/datasets.R), written to bench/proto/r12/review/data/<name>.rds in the
# format of bench/sim.R.   Rscript bench/proto/r12/review/datasets12.R [names...]
#  P_many10  3000 genes x exactly 10 cis SNPs (bench's P_many has ~69: its genes are 3e4 bp apart, cisdist = 1e5)
#  XL_grm_many  N = 2000, GRM, 400 genes x 25 SNPs: per-gene fixed costs at large N (XL_grm has 40 genes)
#  R_real    real-data-like: heavy-tailed SNPs per window (0..~1500, clustered genes, overlapping windows, genes
#            without cis SNPs), GRM, 1-3 samples per individual, 10 covariates, 1% missing genotypes, rare SNPs
source("bench/sim.R")   # ar1_geno(), make_dataset(), DATASETS (main block not run when sourced)
OUT = "bench/proto/r12/review/data"

make_real = function(seed = 7, ni = 400, ng = 200, miss = 0.01, n_cov = 10, effect = 0.4) {
  set.seed(seed)
  inds = sprintf("I%05d", seq_len(ni)); nrep = sample(1:3, ni, TRUE, prob = c(0.2, 0.5, 0.3))
  ind = rep(inds, nrep); N = length(ind)
  cond = unlist(lapply(nrep, function(k) sample(c("ctl", "trt"), k, replace = k > 2)))   # within-individual, uneven
  anno = data.frame(individual = ind, condition = cond, stringsAsFactors = FALSE)
  # genes in clusters: spacing exponential (mean 1e5) with 30% of genes 5-20 kb after the previous one
  gap = ifelse(runif(ng) < 0.3, runif(ng, 5e3, 2e4), rexp(ng, 1 / 1.2e5))
  start = round(1e6 + cumsum(gap)); chr = "chr1"
  geneloc = data.frame(geneid = sprintf("G%06d", seq_len(ng)), chr = chr, left = start, right = start + round(runif(ng, 1e3, 5e4)))
  # SNPs per gene region: heavy tailed (lognormal, median ~40), 8% of regions empty
  spg = pmin(1500L, as.integer(round(exp(rnorm(ng, log(40), 1.1))))); spg[runif(ng) < 0.08] = 0L
  ns = sum(spg); af = runif(ns, 0.005, 0.5)   # includes rare SNPs (filtered by the >= 5 ALT alleles rule)
  G = matrix(0L, ns, ni); pos = numeric(ns); o = 0
  for (g in seq_len(ng)) if (spg[g]) {
    r = o + seq_len(spg[g]); G[r, ] = ar1_geno(spg[g], ni, af[r], runif(1, 0.5, 0.95))
    pos[r] = sort(round(runif(spg[g], geneloc$left[g] - 9e4, geneloc$right[g] + 9e4))); o = o + spg[g]
  }
  storage.mode(G) = "integer"; snp = sprintf("rs%08d", seq_len(ns)); dimnames(G) = list(snp, inds)
  snploc = data.frame(snpid = snp, chr = chr, pos = pos)
  if (miss > 0) G[matrix(runif(length(G)) < miss, nrow(G))] = NA
  trt = as.numeric(cond == "trt")
  covs = matrix(rnorm(N * n_cov), N, dimnames = list(NULL, paste0("pc", seq_len(n_cov))))
  Ui = matrix(rnorm(N * 5), N)
  Y = matrix(rnorm(ng * N), ng) + outer(rnorm(ng, 0, 0.3), trt) + matrix(rnorm(ng * ni, 0, 0.5), ng)[, match(ind, inds)] +
    matrix(rnorm(ng * 5, 0, 0.4), ng) %*% t(Ui) + matrix(rnorm(ng * n_cov, 0, 0.1), ng) %*% t(covs)
  first = cumsum(c(1, spg))[seq_len(ng)]
  causal = sample(which(spg > 0), round(ng * 0.2))
  truth = data.frame(gene = geneloc$geneid[causal], snp = snp[first[causal]], type = ifelse(seq_along(causal) %% 2 == 0, "eQTL", "response"))
  for (k in seq_along(causal)) {
    g = causal[k]; d = G[first[g], ind]; d[is.na(d)] = 0
    Y[g, ] = Y[g, ] + if (truth$type[k] == "eQTL") effect * d else effect * d * trt
  }
  dimnames(Y) = list(geneloc$geneid, paste0("s", seq_len(N)))
  Z = t(scale(t(ar1_geno(3000, ni, runif(3000, 0.05, 0.5), 0)))); Z[is.na(Z)] = 0
  K = crossprod(Z) / nrow(Z); dimnames(K) = list(inds, inds)
  list(input = Y, genotype = G, geneloc = geneloc, snploc = snploc, anno = anno, sample_kernel = K[ind, ind],
       covariates = covs, truth = truth, params = list(type = "real", seed = seed, n_ind = ni, n_gene = ng))
}

D12 = list(
  P_many10    = function() make_dataset(modifyList(base, list(n_ind = 100, n_gene = 3000, snps_per_gene = 10, gene_spacing = 2.5e5, cis_span = 1e4))),
  XL_grm_many = function() make_dataset(modifyList(base, list(n_ind = 1000, n_gene = 400, snps_per_gene = 25, gene_spacing = 2.5e5, cis_span = 2e4, kernel = "grm"))),
  R_real      = function() make_real()
)

if (sys.nframe() == 0) {
  dir.create(OUT, showWarnings = FALSE)
  which = commandArgs(TRUE); if (!length(which)) which = names(D12)
  for (nm in which) {
    d = D12[[nm]](); saveRDS(d, file.path(OUT, paste0(nm, ".rds")), compress = FALSE)
    cat(nm, ": N =", nrow(d$anno), " genes =", nrow(d$geneloc), " SNPs =", nrow(d$snploc), "\n")
  }
}
