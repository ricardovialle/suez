# Out-of-process runner for an installed suez build (a second version of suez cannot be loaded next to the
# baseline): Rscript bench/proto/review/run_lib.R <lib dir> <in.rds> <out.rds>
#   in.rds: list of prepare()-like lists (U, lambda, UX0, C, boot, genes); out.rds: list of run_impl()-like results.
a = commandArgs(TRUE); suppressPackageStartupMessages(library(suez, lib.loc = a[1]))
cases = readRDS(a[2])
saveRDS(lapply(cases, function(P) tryCatch({
  t0 = proc.time()[["elapsed"]]
  out = lapply(P$genes, function(g) { set.seed(g$seed); .Call("suez_fit_gene", P$U, P$lambda, g$Uy, P$UX0, g$G, P$C, P$boot, PACKAGE = "suez") })
  list(out = out, time = proc.time()[["elapsed"]] - t0)
}, error = function(e) e)), a[3])
