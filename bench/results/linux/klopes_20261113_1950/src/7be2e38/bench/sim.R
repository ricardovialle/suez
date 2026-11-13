# Simulated benchmark datasets for map_interaction_qtl (step 2) and learn_covariance (step 1).
# Usage: Rscript bench/sim.R            (writes every dataset of bench/datasets.R to bench/data/)
#        Rscript bench/sim.R <name>...  (only those)
# Each dataset is an RDS list with the arguments of map_interaction_qtl (input, genotype, geneloc,
# snploc, anno, sample_kernel, covariates, local_ancestry) plus `truth` (causal SNPs).
source("bench/datasets.R")

ar1_geno = function(ns, ni, af, rho) {   # SNPs x individuals, LD as an AR(1) latent Gaussian across SNPs
  z = matrix(0, ns, 2 * ni); z[1, ] = rnorm(2 * ni)
  for (j in seq_len(ns)[-1]) z[j, ] = rho * z[j - 1, ] + sqrt(1 - rho^2) * rnorm(2 * ni)
  hap = (pnorm(z) < af) * 1
  hap[, seq(1, 2 * ni, 2)] + hap[, seq(2, 2 * ni, 2)]
}

make_dataset = function(p) {
  set.seed(p$seed)
  ni = p$n_ind; reps = p$reps; ng = p$n_gene; spg = p$snps_per_gene
  inds = sprintf("I%05d", seq_len(ni))
  nrep = if (isTRUE(p$unbalanced)) sample(1:4, ni, TRUE) else rep(reps, ni)   # samples per individual
  ind = rep(inds, nrep); N = length(ind)
  cond = switch(p$condition,
                factor2 = rep(c("ctl", "trt"), length.out = N),
                factor3 = rep(c("a", "b", "c"), length.out = N),
                numeric = rnorm(N),
                none = rep(1, N))
  anno = data.frame(individual = ind, condition = cond, stringsAsFactors = FALSE)
  # genes 1 per 1e5 bp on one chromosome per 2000 genes; SNPs spread over each gene's cis window
  chr = paste0("chr", (seq_len(ng) - 1) %/% 2000 + 1)
  start = ((seq_len(ng) - 1) %% 2000 + 1) * p$gene_spacing
  geneloc = data.frame(geneid = sprintf("G%06d", seq_len(ng)), chr = chr, left = start, right = start + 1000)
  ns = ng * spg
  af = runif(ns, p$af_min, 0.5)
  G = matrix(0L, ns, ni)
  for (g in seq_len(ng)) {   # SNPs of gene g, in LD blocks
    r = (g - 1) * spg + seq_len(spg)
    G[r, ] = ar1_geno(spg, ni, af[r], p$ld_rho)
  }
  storage.mode(G) = "integer"
  snp = sprintf("rs%08d", seq_len(ns))
  dimnames(G) = list(snp, inds)
  snploc = data.frame(snpid = snp, chr = rep(chr, each = spg),
                      pos = rep(start, each = spg) + round(seq(-p$cis_span, p$cis_span, length.out = spg)))
  if (p$miss > 0) G[matrix(runif(length(G)) < p$miss, nrow(G))] = NA
  # expression: condition + individual + latent factors + cis effects
  trt = if (is.numeric(cond)) cond else as.numeric(factor(cond)) - 1
  U = matrix(rnorm(N * p$n_factor), N)
  covs = if (p$n_cov) matrix(rnorm(N * p$n_cov), N, dimnames = list(NULL, paste0("cov", seq_len(p$n_cov)))) else NULL
  Y = matrix(rnorm(ng * N), ng) + outer(rnorm(ng, 0, 0.3), trt) +
    matrix(rnorm(ng * ni, 0, 0.5), ng)[, rep(seq_len(ni), nrep)] +
    matrix(rnorm(ng * p$n_factor, 0, 0.4), ng) %*% t(U)
  causal = sample(ng, round(ng * 0.2))
  truth = data.frame(gene = geneloc$geneid[causal], snp = snp[(causal - 1) * spg + 1],
                     type = ifelse(seq_along(causal) %% 2 == 0, "eQTL", "response"))
  for (k in seq_along(causal)) {
    g = causal[k]; d = G[(g - 1) * spg + 1, ind]; d[is.na(d)] = 0
    Y[g, ] = Y[g, ] + if (truth$type[k] == "eQTL") p$effect * d else p$effect * d * trt
  }
  dimnames(Y) = list(geneloc$geneid, paste0("s", seq_len(N)))
  kernel = switch(p$kernel,
    default = NULL,
    grm = { Z = t(scale(t(ar1_geno(3000, ni, runif(3000, 0.05, 0.5), 0)))); Z[is.na(Z)] = 0
            K = crossprod(Z) / nrow(Z); dimnames(K) = list(inds, inds); K[ind, ind] },
    learned = { L = matrix(rnorm(N * 10), N); same = outer(ind, ind, "==") * 0.3   # low rank + same individual
                tcrossprod(L) / 10 + same },
    lc = { suppressPackageStartupMessages(library(suez, lib.loc = "bench/lib/baseline"))   # step-1 kernel
           learn_covariance(suez:::quantile_normalize(Y), anno, P = 10)$cov_suez })
  list(input = Y, genotype = G, geneloc = geneloc, snploc = snploc, anno = anno, sample_kernel = kernel,
       covariates = covs, truth = truth, params = p)
}

# Tractor-style local-ancestry inputs: 2 ancestries
make_la_dataset = function(p) {
  set.seed(p$seed)
  ni = p$n_ind; reps = p$reps; ng = p$n_gene; spg = p$snps_per_gene; ns = ng * spg
  inds = sprintf("I%05d", seq_len(ni)); ind = rep(inds, each = reps); N = length(ind)
  cond = rep(c("ctl", "trt"), length.out = N)
  theta = rbeta(ni, 4, 2)                       # proportion of ancestry 0
  af = cbind(runif(ns, 0.1, 0.5), runif(ns, 0.1, 0.5))
  D = H = list(matrix(0, ns, ni), matrix(0, ns, ni))
  for (h in 1:2) {
    anc = matrix(runif(ns * ni) > rep(theta, each = ns), ns) + 1   # 1 or 2, SNP x individual (no tract structure)
    for (a in 1:2) {
      on = anc == a
      al = matrix(runif(ns * ni) < af[, a], ns)
      H[[a]] = H[[a]] + on; D[[a]] = D[[a]] + on * al
    }
  }
  snp = sprintf("rs%08d", seq_len(ns))
  for (a in 1:2) dimnames(D[[a]]) = dimnames(H[[a]]) = list(snp, inds)
  names(D) = names(H) = c("anc0", "anc1")
  chr = rep("chr1", ng); start = seq_len(ng) * p$gene_spacing
  geneloc = data.frame(geneid = sprintf("G%06d", seq_len(ng)), chr = chr, left = start, right = start + 1000)
  snploc = data.frame(snpid = snp, chr = "chr1",
                      pos = rep(start, each = spg) + round(seq(-p$cis_span, p$cis_span, length.out = spg)))
  trt = as.numeric(factor(cond)) - 1
  Y = matrix(rnorm(ng * N), ng) + matrix(rnorm(ng * ni, 0, 0.5), ng)[, rep(seq_len(ni), each = reps)]
  for (g in seq_len(ng)[seq(1, ng, 5)]) Y[g, ] = Y[g, ] + 0.5 * D[[1]][(g - 1) * spg + 1, ind] * (1 + trt)
  dimnames(Y) = list(geneloc$geneid, paste0("s", seq_len(N)))
  list(input = Y, genotype = NULL, geneloc = geneloc, snploc = snploc,
       anno = data.frame(individual = ind, condition = cond, stringsAsFactors = FALSE),
       sample_kernel = NULL, covariates = NULL, local_ancestry = list(dosage = D, hapcount = H), params = p)
}

if (sys.nframe() == 0) {
  dir.create("bench/data", showWarnings = FALSE)
  which = commandArgs(TRUE); if (!length(which)) which = names(DATASETS)
  for (nm in which) {
    p = DATASETS[[nm]]
    f = file.path("bench/data", paste0(nm, ".rds"))
    d = if (identical(p$type, "la")) make_la_dataset(p) else make_dataset(p)
    saveRDS(d, f, compress = FALSE)
    cat(nm, ": N =", nrow(d$anno), " genes =", nrow(d$geneloc), " SNPs =", nrow(d$snploc), "\n")
  }
}
