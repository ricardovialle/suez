# Rscript summarize.R <runs.tsv> <base impl>: per dataset x cores, median of paired (same rep) ratios vs base
a = commandArgs(TRUE); r = read.delim(a[1]); b = a[2]
r$key = paste(r$dataset, r$cores, r$rep)
out = NULL
for (im in setdiff(unique(r$impl), b)) {
  p = merge(r[r$impl == b, ], r[r$impl == im, ], by = "key")
  p$cfg = paste(p$dataset.x, p$cores.x)
  s = do.call(rbind, lapply(split(p, p$cfg), function(q) data.frame(impl = im, dataset = q$dataset.x[1], cores = q$cores.x[1], n = nrow(q),
    wall_base = median(q$wall_s.x), wall = median(q$wall_s.y), wall_ratio = round(median(q$wall_s.y / q$wall_s.x), 3),
    rss_base = median(q$maxrss_mb.x), rss = median(q$maxrss_mb.y), rss_ratio = round(median(q$maxrss_mb.y / q$maxrss_mb.x), 3),
    fp_base = median(q$peak_footprint_mb.x), fp = median(q$peak_footprint_mb.y), fp_ratio = round(median(q$peak_footprint_mb.y / q$peak_footprint_mb.x), 3),
    tree_base = median(q$tree_rss_mb.x), tree = median(q$tree_rss_mb.y), tree_ratio = round(median(q$tree_rss_mb.y / q$tree_rss_mb.x), 3))))
  out = rbind(out, s) }
rownames(out) = NULL; print(out[order(out$dataset, out$cores, out$impl), ], row.names = FALSE)
