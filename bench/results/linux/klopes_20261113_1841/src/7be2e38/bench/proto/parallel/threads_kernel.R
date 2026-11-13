# In-process kernel timing of the threaded suez_fit_gene (bench/lib/par_threads) at 1, 2, 4, 8 threads,
# interleaved per rep; checks identical() results across thread counts. Rscript bench/proto/parallel/threads_kernel.R
suppressPackageStartupMessages({ library(suez, lib.loc = "bench/lib/par_threads"); library(Rcpp) })
h = readLines("bench/proto/harness.R"); eval(parse(text = h[!grepl("^suppressPackageStartupMessages", h)]))
fit = function(nt) function(U, lambda, Uy, UX0, G, C, boot) .Call("suez_fit_gene", U, lambda, Uy, UX0, G, C, boot, nt, PACKAGE = "suez")
for (ds in c("S", "M", "M_dense")) {
  P = prepare(ds, n_genes = c(S = 40, M = 30, M_dense = 8)[[ds]])
  tt = matrix(NA, 5, 4, dimnames = list(NULL, paste0("t", c(1, 2, 4, 8)))); ref = NULL; same = TRUE
  for (rep in 1:5) for (nt in c(1, 2, 4, 8)) {
    r = run_impl(fit(nt), P); tt[rep, paste0("t", nt)] = r$time
    if (is.null(ref)) ref = r$out else same = same && identical(ref, r$out)
  }
  med = apply(tt, 2, median); sp = apply(tt[, 1] / tt, 2, median)
  cat(sprintf("%s (%d genes, %d SNPs): median s %s | paired speed-up %s | efficiency %s | identical across threads: %s\n",
              ds, length(P$genes), r$n_snp, paste(round(med, 2), collapse = " / "), paste(round(sp, 2), collapse = " / "),
              paste(round(sp / c(1, 2, 4, 8), 2), collapse = " / "), same))
}
