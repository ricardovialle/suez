# learn_covariance memory: Rscript bench/proto/memory/step1_profile.R <lib> <dataset> [P] [n_genes]
a = commandArgs(TRUE); lib = a[1]; dsn = a[2]; P = if (length(a) > 2) as.integer(a[3]) else 10
suppressPackageStartupMessages({ library(suez, lib.loc = lib); library(doMC) })
rss = function() as.numeric(system(paste("ps -o rss= -p", Sys.getpid()), intern = TRUE)) / 1024
fpk = function() { x = system(paste("footprint -p", Sys.getpid()), intern = TRUE); x[grep("phys_footprint_peak", x)] }
d = readRDS(file.path("bench/data", paste0(dsn, ".rds")))
x = suez:::quantile_normalize(d$input); if (length(a) > 3) x = x[seq_len(min(nrow(x), as.integer(a[4]))), ]
invisible(loadNamespace("rstan")); invisible(gc(reset = TRUE))
cat(sprintf("%s N=%d genes=%d P=%d: before RSS %.0f MB\n", dsn, ncol(x), nrow(x), P, rss()))
set.seed(42); t = system.time(fit <- if (length(a) > 4) learn_covariance(x, d$anno, P = P, iter = as.integer(a[5])) else learn_covariance(x, d$anno, P = P))[["elapsed"]]
g = gc()
cat(sprintf("  after %.1f s: RSS %.0f MB, R heap max %.0f MB, %s\n", t, rss(), g[2, ncol(g)] + g[1, ncol(g)], trimws(fpk())))
