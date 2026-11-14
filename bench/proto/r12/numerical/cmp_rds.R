# compare(a, b) of two saved run_impl results: Rscript cmp_rds.R a.rds b.rds
source("bench/proto/harness.R")
a = commandArgs(TRUE); x = readRDS(a[1]); y = readRDS(a[2]); cm = compare(x, y)
A = do.call(rbind, lapply(x$out, `[[`, "res")); B = do.call(rbind, lapply(y$out, `[[`, "res"))
cat(sprintf("%s: identical %s; dl0 %.1e dl_geno %.1e dl_int %.1e ddf %g dbeta/max(|b|,se) %.1e rel_se %.1e dl_boot %.1e/%.1e NAmis %d\n", basename(a[2]),
  identical(lapply(x$out, unclass), lapply(y$out, unclass)), cm["l0"], cm["l_geno"], cm["l_interact"], cm["df"],
  max(abs(A[, 4] - B[, 4]) / pmax(abs(A[, 4]), A[, 5]), na.rm = TRUE), cm["rel_se"], cm["l_boot_geno"], cm["l_boot_interact"], cm["na_mismatch"]))
