# Interleaved timing + correctness of a prototype vs the installed base12 engine.
#   Rscript bench/proto/r12/numerical/cmp.R <proto.cpp> <reps> <dataset:genes> ...
Sys.setenv(SUEZ_HARNESS_LIB = "bench/lib/base12"); source("bench/proto/harness.R")
a = commandArgs(TRUE); env = new.env(); Rcpp::sourceCpp(a[1], env = env); reps = as.integer(a[2])
for (spec in a[-(1:2)]) {
  ds = sub(":.*", "", spec); ng = as.integer(sub(".*:", "", spec)); P = prepare(ds, ng)
  tb = tp = numeric(reps)
  for (r in seq_len(reps)) { b = run_impl(baseline_fit, P); p = run_impl(env$fit_gene_proto, P); tb[r] = b$time; tp[r] = p$time }
  cm = compare(b, p); bo = do.call(rbind, lapply(b$out, `[[`, "res"))
  relb = max(abs(bo[, 4] - do.call(rbind, lapply(p$out, `[[`, "res"))[, 4]) / pmax(abs(bo[, 4]), bo[, 5]), na.rm = TRUE)
  cat(sprintf("%-8s N=%d snps=%d  base %.3f  proto %.3f  ratio(base/proto) %.3f  [%s]\n", ds, length(P$lambda) + 1, b$n_snp, median(tb), median(tp),
              median(tb / tp), paste(round(tb / tp, 2), collapse = " ")))
  cat(sprintf("   dl0 %.1e dl_geno %.1e dl_int %.1e ddf %g dbeta/max(|b|,se) %.1e rel_se %.1e dl_boot %.1e/%.1e NAmis %d\n", cm["l0"], cm["l_geno"],
              cm["l_interact"], cm["df"], relb, cm["rel_se"], cm["l_boot_geno"], cm["l_boot_interact"], cm["na_mismatch"]))
}
