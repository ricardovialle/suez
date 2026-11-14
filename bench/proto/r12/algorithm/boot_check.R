# bootstrap null LR distribution, base12 vs structured basis (e2e outputs of e2e.sh, rep 1)
for (ds in commandArgs(TRUE)) {
  a = readRDS(sprintf("bench/results/out/algorithm_base_%s_step2_c1_r1.rds", ds)); b = readRDS(sprintf("bench/results/out/algorithm_zz_%s_step2_c1_r1.rds", ds))
  la = 2 * (a$l_boot_interact - a$l_boot_geno); lb = 2 * (b$l_boot_interact - b$l_boot_geno); la = la[is.finite(la)]; lb = lb[is.finite(lb)]
  q = function(x) sprintf("%.3f", quantile(x, c(.5, .9, .99, .999)))
  k = suppressWarnings(ks.test(la, lb))
  cat(sprintf("%s\tn %d/%d  mean %.4f/%.4f  q50,q90,q99,q99.9 %s / %s  KS D %.4f p %.3f\n", ds, length(la), length(lb), mean(la), mean(lb),
    paste(q(la), collapse = ","), paste(q(lb), collapse = ","), k$statistic, k$p.value))
}
