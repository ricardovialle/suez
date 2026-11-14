# Kernel-level check of the threaded fit_gene (pkg_thr) against base12, in one process, interleaved:
# identical() to base12 at every thread count, and time per thread count.
#   Rscript bench/proto/r12/parallel/kernel_threads.R <dataset> <n_genes> <reps> [threads ...]
a = commandArgs(TRUE); dsn = a[1]; ng = as.integer(a[2]); reps = as.integer(a[3])
thr = if (length(a) > 3) as.integer(a[-(1:3)]) else c(1, 2, 4, 8)
Sys.setenv(SUEZ_HARNESS_LIB = "bench/lib/base12"); source("bench/proto/harness.R")
so = file.path(tempdir(), "suez_thr.so")   # renamed: no R_init_suez, so its symbols are looked up dynamically
invisible(file.copy("bench/lib/parallel_thr/suez/libs/suez.so", so, overwrite = TRUE)); dll = dyn.load(so)
sym = getNativeSymbolInfo("suez_fit_gene_threads", dll)
thr_fit = function(nt) function(U, lambda, Uy, UX0, G, C, boot) .Call(sym, U, lambda, Uy, UX0, G, C, boot, as.integer(nt))
P = prepare(dsn, n_genes = ng)
tm = list(); outs = list()
for (r in seq_len(reps)) {
  b = run_impl(baseline_fit, P); tm$base = c(tm$base, b$time); if (r == 1) outs$base = b$out
  for (nt in thr) { k = paste0("t", nt); x = run_impl(thr_fit(nt), P); tm[[k]] = c(tm[[k]], x$time); if (r == 1) outs[[k]] = x$out }
}
same = sapply(outs[-1], function(o) identical(o, outs$base))
med = sapply(tm, median)
cat(sprintf("%s genes=%d snps=%d N=%d | base %.2fs | %s | identical: %s\n", dsn, length(P$genes), b$n_snp, nrow(P$U),
            med[["base"]], paste(sprintf("%s %.2fs (%.2fx)", names(med)[-1], med[-1], med[["base"]] / med[-1]), collapse = ", "),
            paste(names(same), same, collapse = " ")))
