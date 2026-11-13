# Round 11: summarise fit_bench.R outputs. Files are <label>_N<N>_i<0|1>_r<rep>.rds in <dir>.
# Usage: Rscript compare_fit.R <dir> <reference label> [out.tsv]
# Per (N, interaction, label): median elapsed over reps, speed-up vs the reference label, and the largest
# difference of l0 / l_geno / l_interact / beta / se against the reference label, plus NA-pattern mismatches.
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]; ref <- args[2]; out <- if (length(args) > 2) args[3] else file.path(dir, "fit_summary.tsv")
f <- list.files(dir, "^.+_N[0-9]+_i[01]_r[0-9]+\\.rds$")
m <- regmatches(f, regexec("^(.+)_N([0-9]+)_i([01])_r([0-9]+)\\.rds$", f))
key <- data.frame(file = f, label = sapply(m, `[`, 2), N = as.integer(sapply(m, `[`, 3)),
                  inter = as.integer(sapply(m, `[`, 4)), rep = as.integer(sapply(m, `[`, 5)))
rows <- list()
for (g in split(key, list(key$N, key$inter), drop = TRUE)) {
  load1 <- function(k) readRDS(file.path(dir, k$file[1]))
  r0 <- load1(g[g$label == ref, ])
  t_ref <- median(sapply(g$file[g$label == ref], function(x) readRDS(file.path(dir, x))$time[["elapsed"]]))
  for (lab in unique(g$label)) {
    gg <- g[g$label == lab, ]; r <- load1(gg)
    t <- median(sapply(gg$file, function(x) readRDS(file.path(dir, x))$time[["elapsed"]]))
    d <- max(sapply(seq_along(r$res), function(i) max(abs(c(r$res[[i]]$l0 - r0$res[[i]]$l0,
                                                              r$res[[i]]$res[, 1:5] - r0$res[[i]]$res[, 1:5])), na.rm = TRUE)))
    na <- sum(sapply(seq_along(r$res), function(i) sum(is.na(r$res[[i]]$res[, 1:5]) != is.na(r0$res[[i]]$res[, 1:5]))))
    rows[[length(rows) + 1]] <- data.frame(N = gg$N[1], interaction = gg$inter[1], label = lab, version = r$version,
                                           reps = nrow(gg), median_s = t, speedup_vs_ref = t_ref / t,
                                           max_abs_diff_vs_ref = d, na_mismatch = na)
  }
}
res <- do.call(rbind, rows); res <- res[order(res$N, res$interaction, res$label), ]
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
print(res, row.names = FALSE, digits = 3)
