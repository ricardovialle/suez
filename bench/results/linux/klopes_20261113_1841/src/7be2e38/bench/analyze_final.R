# Summary tables and figures of the final benchmark (bench/final_bench.sh -> bench/results/runs_final.tsv).
#   Rscript bench/analyze_final.R  -> bench/results/final_summary.tsv, final_scaling.tsv, *.png
r = read.delim("bench/results/runs_final.tsv")
med = function(x) median(x)
key = c("dataset", "step", "cores")
agg = aggregate(cbind(wall_s, user_s, maxrss_mb, peak_footprint_mb, tree_rss_mb) ~ impl + dataset + step + cores, r, med)
rng = aggregate(wall_s ~ impl + dataset + step + cores, r, function(x) sprintf("%.2f-%.2f", min(x), max(x)))
names(rng)[5] = "wall_range"; agg = merge(agg, rng)
nrep = aggregate(wall_s ~ impl + dataset + step + cores, r, length); names(nrep)[5] = "reps"; agg = merge(agg, nrep)
b = agg[agg$impl == "baseline", ]; f = agg[agg$impl == "final", ]
m = merge(b, f, by = key, suffixes = c("_base", "_final"))
# per-rep paired ratios, to show that the median improvement holds in every repetition
pr = merge(r[r$impl == "baseline", c(key, "rep", "wall_s", "maxrss_mb")], r[r$impl == "final", c(key, "rep", "wall_s", "maxrss_mb")],
           by = c(key, "rep"), suffixes = c("_b", "_f"))
worst = aggregate(cbind(sp = wall_s_b / wall_s_f, ram = maxrss_mb_f / maxrss_mb_b) ~ dataset + step + cores, pr,
                  function(x) c(min = min(x), max = max(x)))
m$speedup = m$wall_s_base / m$wall_s_final
m$runtime_reduction_pct = 100 * (1 - m$wall_s_final / m$wall_s_base)
m$cpu_speedup = m$user_s_base / m$user_s_final
m$rss_reduction_pct = 100 * (1 - m$maxrss_mb_final / m$maxrss_mb_base)
m$tree_rss_reduction_pct = 100 * (1 - m$tree_rss_mb_final / m$tree_rss_mb_base)
m = merge(m, data.frame(worst[key], speedup_min = worst$sp[, "min"], speedup_max = worst$sp[, "max"],
                        ram_ratio_max = worst$ram[, "max"]), by = key)
m = m[order(m$step, m$cores, -m$wall_s_base), ]
out = m[, c(key, "reps_base", "wall_s_base", "wall_range_base", "wall_s_final", "wall_range_final", "speedup", "speedup_min",
            "speedup_max", "runtime_reduction_pct", "cpu_speedup", "maxrss_mb_base", "maxrss_mb_final", "rss_reduction_pct",
            "tree_rss_mb_base", "tree_rss_mb_final", "tree_rss_reduction_pct")]
write.table(out, "bench/results/final_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 250); print(format(out, digits = 3), row.names = FALSE)

# scaling with cores (M, M_dense)
sc = agg[agg$step == "step2" & agg$dataset %in% c("M", "M_dense"), ]
sc$speedup_vs_1core = 1 / (sc$wall_s / sc$wall_s[match(paste(sc$impl, sc$dataset, 1), paste(sc$impl, sc$dataset, sc$cores))])
sc$efficiency = sc$speedup_vs_1core / sc$cores
write.table(sc[order(sc$impl, sc$dataset, sc$cores), c("impl", "dataset", "cores", "wall_s", "tree_rss_mb", "speedup_vs_1core", "efficiency")],
            "bench/results/final_scaling.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# figures
one = m[m$step == "step2" & m$cores == 1, ]; one = one[order(one$wall_s_base), ]
png("bench/results/final_runtime.png", 1500, 900, res = 150)
par(mar = c(7, 5, 3, 1))
bp = barplot(rbind(one$wall_s_base, one$wall_s_final), beside = TRUE, names.arg = one$dataset, las = 2, log = "y",
             col = c("grey65", "steelblue"), ylab = "wall-clock (s, log scale)", main = "Step 2 runtime, 1 core (median of 3)",
             ylim = c(0.4, max(one$wall_s_base) * 2.2))
text(colMeans(bp), pmax(one$wall_s_base, one$wall_s_final) * 1.25, sprintf("%.1fx", one$speedup), cex = 0.7)
legend("topleft", c("0.4 (baseline)", "optimized"), fill = c("grey65", "steelblue"), bty = "n")
dev.off()
png("bench/results/final_ram.png", 1500, 900, res = 150)
par(mar = c(7, 5, 3, 1))
barplot(rbind(one$maxrss_mb_base, one$maxrss_mb_final), beside = TRUE, names.arg = one$dataset, las = 2,
        col = c("grey65", "steelblue"), ylab = "peak RSS (MB)", main = "Step 2 peak RAM, 1 core (median of 3)")
legend("topleft", c("0.4 (baseline)", "optimized"), fill = c("grey65", "steelblue"), bty = "n")
dev.off()
png("bench/results/final_scaling.png", 1500, 700, res = 150)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))
for (ds in c("M", "M_dense")) {
  s = sc[sc$dataset == ds, ]
  plot(NA, xlim = c(1, 8), ylim = c(0, max(s$wall_s)), xlab = "cores", ylab = "wall-clock (s)", main = ds, xaxt = "n")
  axis(1, c(1, 2, 4, 8))
  for (i in c("baseline", "final")) { z = s[s$impl == i, ]; z = z[order(z$cores), ]
    lines(z$cores, z$wall_s, type = "b", pch = 19, col = if (i == "final") "steelblue" else "grey40") }
  legend("topright", c("0.4", "optimized"), col = c("grey40", "steelblue"), lwd = 2, bty = "n")
}
dev.off()
png("bench/results/final_tree_ram.png", 1500, 700, res = 150)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))
for (ds in c("M", "M_dense")) {
  s = sc[sc$dataset == ds, ]
  plot(NA, xlim = c(1, 8), ylim = c(0, max(s$tree_rss_mb)), xlab = "cores", ylab = "peak RSS, process tree (MB)", main = ds, xaxt = "n")
  axis(1, c(1, 2, 4, 8))
  for (i in c("baseline", "final")) { z = s[s$impl == i, ]; z = z[order(z$cores), ]
    lines(z$cores, z$tree_rss_mb, type = "b", pch = 19, col = if (i == "final") "steelblue" else "grey40") }
  legend("topleft", c("0.4", "optimized"), col = c("grey40", "steelblue"), lwd = 2, bty = "n")
}
dev.off()
