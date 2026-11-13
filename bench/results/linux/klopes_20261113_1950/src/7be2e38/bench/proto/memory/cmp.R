# Rscript bench/proto/memory/cmp.R a.rds b.rds : identical? else max abs diff per numeric column
a = readRDS(commandArgs(TRUE)[1]); b = readRDS(commandArgs(TRUE)[2])
if (identical(a, b)) cat("IDENTICAL\n") else {
  cat("not identical; dims", dim(a), "/", dim(b), " names equal", identical(names(a), names(b)), " class", class(a), "/", class(b), "\n")
  for (n in intersect(names(a), names(b))) if (is.numeric(a[[n]])) cat(n, suppressWarnings(max(abs(a[[n]] - b[[n]]), na.rm = TRUE)), sum(is.na(a[[n]]) != is.na(b[[n]])), "\n")
    else if (!identical(a[[n]], b[[n]])) cat(n, "differs\n")
}
