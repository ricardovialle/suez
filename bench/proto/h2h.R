# Head-to-head kernel-level timing and correctness: p0 vs candidate prototypes (interleaved).
# Usage: Rscript bench/proto/h2h.R <n_genes> <dataset>... ; candidates listed in CANDS below.
source("bench/proto/harness.R"); source("bench/proto/review/interleave.R")
CANDS = c(p0 = "bench/proto/p0_baseline.cpp", p3 = "bench/proto/algorithm/p3_compress_fast.cpp",
          p7 = "bench/proto/numerical/p7_exact.cpp")
a = commandArgs(TRUE); ng = as.integer(a[1]); dss = a[-1]
f = lapply(CANDS, function(file) { e = new.env(); Rcpp::sourceCpp(file, env = e); e$fit_gene_proto })
rel_se_beta = function(a, b) {   # |dbeta| / max(|beta|, se), the reviewer's criterion (H8)
  A = do.call(rbind, lapply(a$out, `[[`, "res")); B = do.call(rbind, lapply(b$out, `[[`, "res"))
  max(abs(A[, 4] - B[, 4]) / pmax(abs(A[, 4]), A[, 5]), na.rm = TRUE) }
rows = list()
for (ds in dss) {
  P = prepare(ds, ng)
  base = run_impl(baseline_fit, P)
  t = time_interleaved(f, P, reps = 3)$summary
  for (nm in names(f)[-1]) {
    cm = compare(base, run_impl(f[[nm]], P))
    rows[[length(rows) + 1]] = data.frame(dataset = ds, impl = nm, n = length(P$lambda),
      cpu_p0 = t$median_cpu[1], cpu = t$median_cpu[t$impl == nm], speedup = 1 / t$ratio_cpu[t$impl == nm],
      max_dl = max(cm[c("l0", "l_geno", "l_interact")]), max_dl_boot = max(cm[c("l_boot_geno", "l_boot_interact")]),
      dbeta_se = rel_se_beta(base, run_impl(f[[nm]], P)), rel_se = cm[["rel_se"]], na_mismatch = cm[["na_mismatch"]])
  }
}
r = do.call(rbind, rows); print(r, digits = 3, row.names = FALSE)
write.table(r, "bench/results/h2h_round1.tsv", sep = "\t", quote = FALSE, row.names = FALSE,
            append = file.exists("bench/results/h2h_round1.tsv"), col.names = !file.exists("bench/results/h2h_round1.tsv"))
