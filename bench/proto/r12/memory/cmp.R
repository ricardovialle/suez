# Rscript cmp.R <implA> <implB> [suffix]: identical() of every out/<A>_*<suffix> with out/<B>_*
a = commandArgs(TRUE); suf = if (length(a) >= 3) a[3] else "_ck.rds"
fs = list.files("out", pattern = paste0("^", a[1], "_.*", suf, "$"))
for (f in fs) { g = file.path("out", sub(paste0("^", a[1]), a[2], f)); if (!file.exists(g)) next
  x = readRDS(file.path("out", f)); y = readRDS(g); num = sapply(x, is.numeric)
  md = if (identical(dim(x), dim(y))) max(sapply(names(x)[num], function(k) max(abs(x[[k]] - y[[k]]), 0, na.rm = TRUE))) else NA
  cat(sprintf("%-40s identical=%s rows=%d maxdiff=%s\n", f, identical(x, y), nrow(x), md)) }
