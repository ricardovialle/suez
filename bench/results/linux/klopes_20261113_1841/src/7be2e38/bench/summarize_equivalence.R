# Data-column equivalence of a build vs the baseline on every dataset (bootstrap columns excluded: their draws
# differ by design since round 7). Rscript bench/summarize_equivalence.R <impl> -> bench/results/equivalence_<impl>.tsv
impl = commandArgs(TRUE)[1]
DS = c("S", "M", "L", "XL", "M_grm", "M_learn", "M_num", "M_fac3", "M_eqtl", "M_cov", "M_miss", "M_dense", "S_rare",
       "LA", "L_grm", "M_unbal", "M_lc", "HS")
rows = lapply(DS, function(d) {
  a = as.data.frame(readRDS(sprintf("bench/results/out/baseline_%s_step2_c1_r0.rds", d)))
  b = as.data.frame(readRDS(sprintf("bench/results/out/%s_%s_step2_c1_r0.rds", impl, d)))
  stopifnot(identical(paste(a$gene, a$cis_snp), paste(b$gene, b$cis_snp)))
  ll = setdiff(grep("^l0$|^l_", names(a), value = TRUE), grep("boot", names(a), value = TRUE))
  pc = grep("^p_", setdiff(names(a), grep("boot", names(a), value = TRUE)), value = TRUE)
  be = grep("^beta", names(a), value = TRUE); se = sub("^beta", "se", be)
  num = c(ll, pc, be, se, grep("^df", names(a), value = TRUE))
  na_mis = sum(is.na(as.matrix(a[num])) != is.na(as.matrix(b[num])))
  dl = max(abs(as.matrix(a[ll]) - as.matrix(b[ll])), na.rm = TRUE)
  lp = function(x) -log10(pmax(x, 1e-300))
  dp = max(abs(lp(as.matrix(a[pc])) - lp(as.matrix(b[pc]))), na.rm = TRUE)
  db = max(abs(as.matrix(a[be]) - as.matrix(b[be])) / pmax(abs(as.matrix(a[be])), as.matrix(a[se])), na.rm = TRUE)
  ds = max(abs(as.matrix(b[se]) / as.matrix(a[se]) - 1), na.rm = TRUE)
  lead = function(x, p) tapply(seq_len(nrow(x)), x$gene, function(i) x$cis_snp[i][which.min(x[[p]][i])])
  same_lead = all(sapply(intersect(pc, c("p_geno", "p_interact")), function(p) identical(lead(a, p), lead(b, p))))
  sig = function(x, t) sum(as.matrix(x[pc]) < t, na.rm = TRUE)
  data.frame(dataset = d, tests = nrow(a), max_abs_dl = dl, max_dlog10p = dp, max_dbeta_over_se = db, max_rel_dse = ds,
             na_mismatch = na_mis, lead_snps_identical = same_lead, sig_1e5_base = sig(a, 1e-5), sig_1e5_new = sig(b, 1e-5),
             pass = dl <= 1e-9 && db <= 1e-6 && ds <= 1e-6 && na_mis == 0 && same_lead)
})
r = do.call(rbind, rows); print(r, digits = 3, row.names = FALSE)
write.table(r, sprintf("bench/results/equivalence_%s.tsv", impl), sep = "\t", quote = FALSE, row.names = FALSE)
