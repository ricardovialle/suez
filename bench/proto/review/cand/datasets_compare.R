# Kernel-level comparison of candidates vs baseline on benchmark datasets (first n genes), next to the
# baseline's own rounding-noise floor (baseline on Uy * (1 + 1e-15 N(0,1)), as bench/proto/algorithm/floor.R).
#   Rscript bench/proto/review/cand/datasets_compare.R cand1.cpp [cand2.cpp ...]
source("bench/proto/harness.R"); source("bench/proto/review/validate.R")
files = commandArgs(TRUE); impl = list()
for (f in files) { e = new.env(); Rcpp::sourceCpp(f, env = e); impl[[sub("\\.cpp$", "", basename(f))]] = e$fit_gene_proto }
dss = c(S = 30, M_grm = 20, M_learn = 20, M_lc = 20, M_fac3 = 20, M_eqtl = 20, M_cov = 20, M_unbal = 20, L_grm = 8, HS = 30, S_rare = 30, M_num = 20)
summ = function(a, b, P) {
  A = harness_table(a, P); B = harness_table(b, P); v = validate(A, B, quiet = TRUE)
  cl = v$columns; mx = function(cs) suppressWarnings(max(cl$max_abs[cl$column %in% cs], na.rm = TRUE))
  db = max(abs(A$beta_geno - B$beta_geno) / A$se_geno, na.rm = TRUE)
  data.frame(dl_data = mx(c("l0", "l_geno", "l_interact")), dl_boot = mx(c("l_boot_geno", "l_boot_interact")),
             dbeta_se = db, na_mis = sum(cl$na_mismatch), lead_same = min(v$ranking$lead_same / v$ranking$genes))
}
out = list()
for (ds in names(dss)) {
  P = prepare(ds, dss[[ds]]); base = run_impl(baseline_fit, P)
  Pf = P; set.seed(3); for (i in seq_along(Pf$genes)) Pf$genes[[i]]$Uy = Pf$genes[[i]]$Uy * (1 + 1e-15 * rnorm(length(Pf$genes[[i]]$Uy)))
  fl = run_impl(baseline_fit, Pf)
  rows = list(cbind(ds = ds, impl = "floor(baseline 1e-15)", summ(base, fl, P), t_ratio = NA))
  for (k in names(impl)) { r = run_impl(impl[[k]], P); rows[[k]] = cbind(ds = ds, impl = k, summ(base, r, P), t_ratio = base$time / r$time) }
  out[[ds]] = do.call(rbind, rows); print(format(out[[ds]], digits = 3), row.names = FALSE)
}
saveRDS(do.call(rbind, out), "bench/proto/review/out/cand_datasets_compare.rds")
