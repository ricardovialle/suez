# Rscript bind_mem.R <mode: none|bind_rows|lean> [genes] [snps]: peak RSS of combining per-gene result tables
GCK = nzchar(Sys.getenv("GCK")); a = commandArgs(TRUE); mode = a[1]; ng = as.integer(a[2]); S = as.integer(a[3])
snps = sprintf("chr1:%d:A:G", 1:(ng * S / 5))   # SNPs shared by ~5 genes, as overlapping cis windows
mk = function(i) { m = matrix(rnorm(S * 7), S); gene = sprintf("ENSG%011d", i); l0 = -100
  structure(list(gene = rep(gene, S), cis_snp = snps[((i - 1) * S / 5 + seq_len(S) - 1) %% length(snps) + 1], l0 = rep(l0, S),
    l_geno = m[, 1], l_interact = m[, 2], p_geno = pchisq(pmax(0, 2 * (m[, 1] - l0)), 1, lower.tail = FALSE),
    df_interact = rep(1L, S), p_interact = m[, 3], beta_geno = m[, 4], se_geno = m[, 5], l_boot_geno = m[, 6],
    l_boot_interact = m[, 7]), class = "data.frame", row.names = c(NA, -S)) }
set.seed(1); res = lapply(seq_len(ng), mk); invisible(gc())
rss = function() as.numeric(system(sprintf("ps -o rss= -p %d", Sys.getpid()), intern = TRUE)) / 1024
r0 = rss()
# lean: one column at a time, dropping it from the per-gene tables (peak ~ one copy + one column)
bind_lean = function(res) {
  res = res[!vapply(res, is.null, NA)]; nm = names(res[[1]]); n = sum(vapply(res, function(x) length(x[[1]]), 0L))
  out = vector("list", length(nm)); names(out) = nm
  for (k in nm) { out[[k]] = unlist(lapply(res, `[[`, k), use.names = FALSE); res = lapply(res, function(x) { x[[k]] = NULL; x }); if (GCK) gc(FALSE) }
  structure(out, class = "data.frame", row.names = c(NA, -n)) }
t = system.time(out <- switch(mode, none = NULL, bind_rows = dplyr::bind_rows(res), lean = bind_lean(res)))[[3]]
if (mode == "check") stopifnot(identical(dplyr::bind_rows(res), bind_lean(res)))
cat(sprintf("%s genes=%d rows=%d res_list=%.0f MB rss_before=%.0f MB time=%.2f s\n", mode, ng, ng * S, object.size(res) / 2^20, r0, t))
