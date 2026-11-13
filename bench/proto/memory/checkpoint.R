# checkpoint_dir cost: Rscript bench/proto/memory/checkpoint.R <lib> <dataset> <dir>
a = commandArgs(TRUE); suppressPackageStartupMessages({ library(suez, lib.loc = a[1]); library(doMC) }); registerDoMC(1)
d = readRDS(file.path("bench/data", paste0(a[2], ".rds"))); dir = a[3]; unlink(dir, recursive = TRUE)
run = function(cp) { set.seed(42); system.time(r <- map_interaction_qtl(d$input, d$genotype, d$geneloc, d$snploc, d$anno,
  sample_kernel = d$sample_kernel, checkpoint_dir = cp))[["elapsed"]] -> t; list(r = r, t = t) }
r0 = run(NULL); r1 = run(dir); r2 = run(dir)
f = list.files(dir, full.names = TRUE)
cat(sprintf("%s: no checkpoint %.2f s; writing %.2f s; resuming from %d files (%.1f MB) %.2f s\n", a[2], r0$t, r1$t, length(f), sum(file.size(f)) / 2^20, r2$t))
cat("write run identical to no-checkpoint run:", identical(r0$r, r1$r), "\n")
cat("resumed run: class", class(r2$r), "; max |diff| numeric cols:",
  max(sapply(names(r0$r)[sapply(r0$r, is.numeric)], function(n) max(abs(r0$r[[n]] - r2$r[[n]]), na.rm = TRUE))), "; column types differ:",
  paste(names(r0$r)[sapply(names(r0$r), function(n) !identical(class(r0$r[[n]]), class(r2$r[[n]])))], collapse = ","), "\n")
