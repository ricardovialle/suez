# Medians of bench/proto/r12/parallel/runs.tsv per configuration, speed-up vs the 1-core baseline, parallel
# efficiency (speed-up / CPUs = forks x threads), and identical() of every output vs the 1-core baseline's.
#   Rscript bench/proto/r12/parallel/summarize.R
r = read.delim("bench/proto/r12/parallel/runs.tsv", stringsAsFactors = FALSE)
r$threads = as.integer(ifelse(grepl("_t[0-9]+", r$impl), sub(".*_t([0-9]+).*", "\\1", r$impl), "1"))
r$cpu = r$user_s + r$sys_s
s = aggregate(cbind(wall_s, cpu, tree_rss_mb, maxrss_mb) ~ impl + dataset + cores + threads, r, median)
s$n = aggregate(wall_s ~ impl + dataset + cores + threads, r, length)$wall_s
b = s[s$impl == "parallel_base12" & s$cores == 1, c("dataset", "wall_s")]; names(b)[2] = "base1"
s = merge(s, b, all.x = TRUE)
s$speedup = round(s$base1 / s$wall_s, 2); s$cpus = s$cores * s$threads; s$eff = round(s$speedup / s$cpus, 2)
out = function(impl, ds, cores) file.path("bench/results/out", sprintf("%s_%s_step2_c%d_rr1.rds", impl, ds, cores))
s$identical = mapply(function(i, d, c) {
  a = out("parallel_base12", d, 1L); x = out(i, d, c)
  if (!file.exists(a) || !file.exists(x)) NA else identical(readRDS(a), readRDS(x))
}, s$impl, s$dataset, s$cores)
s = s[order(s$dataset, s$wall_s), c("dataset", "impl", "cores", "threads", "cpus", "n", "wall_s", "cpu", "tree_rss_mb",
                                     "maxrss_mb", "speedup", "eff", "identical")]
print(s, row.names = FALSE)
write.table(s, "bench/proto/r12/parallel/summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
