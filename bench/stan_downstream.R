# Step-1 Stan rewrite (bench/proto/memory/stan.diff): effect of its kernel on step-2 results.
# Rscript bench/stan_downstream.R <dataset>
a = commandArgs(TRUE); ds = a[1]
d = readRDS(file.path("bench/data", paste0(ds, ".rds")))
kern = list()
for (lib in c("baseline", "mem_stan")) {
  f = tempfile(fileext = ".rds")
  system2("Rscript", c("-e", shQuote(sprintf(
    "suppressPackageStartupMessages(library(suez, lib.loc='bench/lib/%s')); library(doMC); registerDoMC(1); d=readRDS('bench/data/%s.rds'); set.seed(42); t0=proc.time()[[3]]; lc=learn_covariance(suez:::quantile_normalize(d$input), d$anno, P=10); t1=proc.time()[[3]]-t0; g=head(d$geneloc$geneid,40); set.seed(1); r=map_interaction_qtl(d$input[g,], d$genotype, d$geneloc, d$snploc, d$anno, sample_kernel=lc$cov_suez, permutation_approach='none'); saveRDS(list(K=lc$cov_suez, lp=lc$optimized_fit$value, t=t1, r=r), '%s')", lib, ds, f))), stdout = FALSE, stderr = FALSE)
  kern[[lib]] = readRDS(f)
}
b = kern$baseline; s = kern$mem_stan
cat(ds, "step1 time baseline", b$t, "s, rewrite", s$t, "s; log density", b$lp, "vs", s$lp, "\n")
cat("cov_suez max rel diff", max(abs(b$K - s$K)) / max(abs(b$K)), "\n")
lp = function(p) -log10(p)
cat("max |d -log10 p_geno|", max(abs(lp(b$r$p_geno) - lp(s$r$p_geno))), " p_interact", max(abs(lp(b$r$p_interact) - lp(s$r$p_interact))), "\n")
lead = function(r) tapply(seq_len(nrow(r)), r$gene, function(i) r$cis_snp[i][which.min(r$p_geno[i])])
cat("lead SNP identical in", mean(lead(b$r) == lead(s$r)), "of genes; Spearman p_geno", cor(b$r$p_geno, s$r$p_geno, method = "spearman"), "\n")
