# Medians (and min-max) of my runs in bench/results/runs.tsv, speed-up vs baseline 1 core, and identical()
# of every prototype output against the baseline 1-core output.
r = read.delim("bench/results/runs.tsv", stringsAsFactors = FALSE)
r = r[grepl("^(p|v)[0-9]$", r$rep), ]
r$label = ifelse(grepl("^par_thr", r$impl), sprintf("%s (cores %d x threads %s)", r$impl, r$cores, sub("par_thr", "", r$impl)), sprintf("%s (cores %d)", r$impl, r$cores))
base1 = tapply(r$wall_s[r$impl == "baseline" & r$cores == 1], r$dataset[r$impl == "baseline" & r$cores == 1], median)
agg = do.call(rbind, lapply(split(r, list(r$dataset, r$label), drop = TRUE), function(x) {
  units = x$cores[1] * if (grepl("^par_thr", x$impl[1])) as.integer(sub("par_thr", "", x$impl[1])) else 1
  w = median(x$wall_s)
  data.frame(dataset = x$dataset[1], variant = x$label[1], cpus = units, n = nrow(x), wall_med = w,
             wall_range = sprintf("%.2f-%.2f", min(x$wall_s), max(x$wall_s)), user_med = median(x$user_s),
             sys_med = median(x$sys_s), util = median(x$cpu_util), speedup = base1[x$dataset[1]] / w,
             efficiency = base1[x$dataset[1]] / w / units)
}))
agg = agg[order(agg$dataset, agg$cpus, agg$variant), ]
print(format(agg, digits = 3), row.names = FALSE)
# within-rep paired ratio vs baseline at the same cores (reps v1..v3)
v = r[grepl("^v", r$rep), ]
b = v[v$impl == "baseline", c("dataset", "cores", "rep", "wall_s")]; names(b)[4] = "b"
m = merge(v[v$impl != "baseline" & !grepl("^par_thr", v$impl), ], b)
if (nrow(m)) { cat("\npaired ratio baseline/variant (same rep, same cores), median [min, max]:\n")
  print(aggregate(cbind(ratio = b / wall_s) ~ dataset + cores + impl, m, function(z) sprintf("%.2f [%.2f, %.2f]", median(z), min(z), max(z)))) }
o = "bench/results/out/"; f = list.files(o, "^(par_|baseline_).*_r(p|v|chk)[0-9a-z]*\\.rds$")
ok = sapply(f, function(x) { ds = sub("^(par_[a-z0-9]+|baseline)_(.*)_step2_.*", "\\2", x)
  ref = paste0(o, "baseline_", ds, "_step2_c1_rp1.rds"); if (!file.exists(ref)) NA else identical(readRDS(ref), readRDS(paste0(o, x))) })
cat("\nidentical() to baseline 1-core output:", sum(ok, na.rm = TRUE), "of", sum(!is.na(ok)), "files;", "not identical:", paste(names(ok)[which(!ok)], collapse = ", "), "\n")
