# Distributional comparison of bootstrap null LRs between two builds (used when a change alters which draws
# a seed produces, e.g. the eigen basis): Rscript bench/boot_distribution.R <ref.rds> <cand.rds>
a = commandArgs(TRUE); A = as.data.frame(readRDS(a[1])); B = as.data.frame(readRDS(a[2]))
lr = function(d) if ("l_boot_interact_null" %in% names(d)) 2 * (d$l_boot_interact - d$l_boot_interact_null) else 2 * (d$l_boot_interact - d$l_boot_geno)
x = lr(A); y = lr(B); x = x[is.finite(x)]; y = y[is.finite(y)]
q = c(0.5, 0.9, 0.99, 0.999)
ks = suppressWarnings(ks.test(x, y))
pb = function(d) { l = lr(d); o = 2 * (d$l_interact - if ("l_interact_null" %in% names(d)) d$l_interact_null else d$l_geno)
  b = sort(l[is.finite(l)]); 1 - findInterval(o, b) / length(b) }   # pooled bootstrap p per test
pa = pb(A); pbb = pb(B); ok = is.finite(pa) & is.finite(pbb)
cat(sprintf("n %d/%d  mean %.4f/%.4f  q50 %.4f/%.4f  q90 %.3f/%.3f  q99 %.3f/%.3f  q99.9 %.3f/%.3f  KS D %.4f p %.3g  pooled-boot p: max |d -log10 p| %.3f, n(p<1e-3) %d/%d\n",
  length(x), length(y), mean(x), mean(y), quantile(x, .5), quantile(y, .5), quantile(x, .9), quantile(y, .9),
  quantile(x, .99), quantile(y, .99), quantile(x, .999), quantile(y, .999), ks$statistic, ks$p.value,
  max(abs(log10(pa[ok] + 1e-12) - log10(pbb[ok] + 1e-12))), sum(pa < 1e-3, na.rm = TRUE), sum(pbb < 1e-3, na.rm = TRUE)))
