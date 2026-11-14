# Round-12 C++ phase profile of the baseline engine (p12_timed.cpp) on per-gene inputs, and a check that the
# sourceCpp copy (p12_base.cpp) reproduces the installed baseline exactly.
#   Rscript bench/proto/r12/profile12.R [datasets...]      (from the repo root)
Sys.setenv(SUEZ_HARNESS_LIB = "bench/lib/base12")
source("bench/proto/harness.R")
Rcpp::sourceCpp("bench/proto/r12/p12_timed.cpp")
a = commandArgs(TRUE)
DS = if (length(a)) a else c("S", "M", "L", "XL", "XXL", "M_grm", "L_grm", "XL_grm", "M_fac3", "M_eqtl", "M_cov", "M_cov30", "M_dense", "P_many")
NG = c(S = 40, M = 30, L = 15, XL = 8, XXL = 4, M_grm = 30, L_grm = 12, XL_grm = 6, M_fac3 = 30, M_eqtl = 30,
       M_cov = 30, M_cov30 = 20, M_dense = 10, P_many = 300)
rows = list()
for (ds in DS) {
  P = prepare(ds, NG[[ds]])
  b = run_impl(baseline_fit, P); timers(TRUE)
  r = run_impl(fit_gene_proto, P); tt = timers(TRUE)
  cmp = compare(b, r)
  sh = round(100 * tt / r$time, 1)
  cat(sprintf("%-8s N=%4d snps=%5d  C++ %.2fs (installed %.2fs)  identical=%s\n", ds, length(P$lambda) + 1, r$n_snp, r$time, b$time,
              identical(lapply(b$out, unclass), lapply(r$out, unclass))))
  print(sh)
  rows[[ds]] = data.frame(dataset = ds, N = length(P$lambda) + 1, snps = r$n_snp, genes = length(P$genes), total_s = r$time, t(tt))
}
out = do.call(rbind, rows)
write.table(out, "bench/results/r12_profile_base.tsv", sep = "\t", quote = FALSE, row.names = FALSE, append = file.exists("bench/results/r12_profile_base.tsv"),
            col.names = !file.exists("bench/results/r12_profile_base.tsv"))
