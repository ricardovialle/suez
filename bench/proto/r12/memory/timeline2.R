# Rscript timeline2.R <lib> <dataset> [n_genes]: phys_footprint now / peak so far (footprint(1)) at each step of map_interaction_qtl
a = commandArgs(TRUE); suppressPackageStartupMessages({ library(suez, lib.loc = a[1]); library(doMC) }); registerDoMC(1)
fp = function(s) { o = system(sprintf("footprint -p %d --noCategories -f bytes", Sys.getpid()), intern = TRUE)
  v = as.numeric(sub(".*: ([0-9]+) B.*", "\\1", grep("phys_footprint", o, value = TRUE))) / 2^20
  cat(sprintf("%-34s footprint %7.1f MB   peak so far %7.1f MB\n", s, v[1], v[2])) }
fp("loaded")
d = readRDS(file.path("/Users/rvialle/suez/bench/data", paste0(a[2], ".rds"))); fp("data read")
ns = asNamespace("suez")
trace("kernel_contrast_eigen", quote(fp("enter kernel_contrast_eigen")), exit = quote(fp("exit kernel_contrast_eigen")), print = FALSE, where = ns)
trace("check_identifiable", quote(fp("set-up done (check_identifiable)")), print = FALSE, where = ns)
trace("pchisq", exit = quote(if (!is.null(df) && length(df) > 1) fp("gene fitted (pchisq)")), print = FALSE, where = ns)
ng = if (length(a) >= 3) as.integer(a[3]) else 4; g = head(d$geneloc$geneid, ng)
set.seed(42); r = map_interaction_qtl(d$input[g, , drop = FALSE], d$genotype, d$geneloc, d$snploc, d$anno, sample_kernel = d$sample_kernel)
fp(sprintf("done, %d genes", ng)); rm(r); invisible(gc()); fp("after gc")
