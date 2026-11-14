# Edge-case suite of bench/proto/review/edge_cases.R (its 23+ inputs, unchanged) for an installed candidate build,
# compared strictly with the reference build: NA and +-Inf patterns identical (edge_cases.R let an Inf-vs-finite
# difference pass), l <= 1e-9 PASS, <= 1e-6 WARN (the round-1 tolerance for extreme inputs), bootstrap <= 1e-5,
# beta/se relative <= 1e-6, df identical, both-error = PASS.
#   Rscript bench/proto/r12/review/edge12.R <cand lib dir> <ref lib dir> <out.tsv>
a = commandArgs(TRUE); cand_lib = a[1]; Sys.setenv(SUEZ_HARNESS_LIB = a[2]); out = a[3]
src = readLines("bench/proto/review/edge_cases.R")
src = src[seq_len(grep("^options\\(width = 200\\)", src)[1] - 1)]   # case construction only
src = sub("a = commandArgs(TRUE)", "a = 'lib:'", src, fixed = TRUE)   # no sourceCpp: candidates are installed builds
eval(parse(text = src))                                               # -> cases, tag, harness (reference lib)
fi = tempfile(fileext = ".rds"); fo = tempfile(fileext = ".rds"); saveRDS(cases, fi)
system2("Rscript", c("bench/proto/review/run_lib.R", cand_lib, fi, fo)); ext = readRDS(fo)
rows = lapply(names(cases), function(nm) {
  b = tryCatch(run_impl(baseline_fit, cases[[nm]]), error = function(e) e); c = ext[[nm]]
  row = function(v, st, det) data.frame(dataset = nm, check = paste("edge case", tag[[nm]]), value = v, threshold = "l 1e-9 (WARN 1e-6), rel 1e-6",
                                        status = st, detail = det)
  if (inherits(b, "error") || inherits(c, "error")) {
    both = inherits(b, "error") && inherits(c, "error")
    return(row(NA, if (both) "PASS" else "FAIL", if (both) "both error" else paste("error:", conditionMessage(if (inherits(c, "error")) c else b))))
  }
  A = do.call(rbind, lapply(b$out, `[[`, "res")); B = do.call(rbind, lapply(c$out, `[[`, "res"))
  A = cbind(sapply(b$out, `[[`, "l0")[rep(seq_along(b$out), sapply(b$out, function(o) nrow(o$res)))], A)
  B = cbind(sapply(c$out, `[[`, "l0")[rep(seq_along(c$out), sapply(c$out, function(o) nrow(o$res)))], B)
  # columns: l0 l_geno l_interact df beta se l_boot_geno l_boot_interact
  pat = sum(is.na(A) != is.na(B) | is.infinite(A) != is.infinite(B))
  ok = is.finite(A) & is.finite(B); D = abs(A - B); D[!ok] = 0
  dl = max(D[, 1:3]); db = max(D[, 7:8])
  rel = function(j) { m = pmax(abs(A[, j]), abs(B[, j]), if (j == 5) A[, 6] else 0, na.rm = TRUE); max(c(0, ifelse(D[, j] == 0, 0, D[, j] / m)), na.rm = TRUE) }
  df_same = identical(A[, 4], B[, 4])
  hard = pat > 0 || !df_same || rel(5) > 1e-6 || rel(6) > 1e-6 || db > 1e-5 || dl > 1e-6
  st = if (hard) (if (tag[[nm]] == "report") "WARN" else "FAIL") else if (dl > 1e-9) "WARN" else "PASS"
  row(dl, st, sprintf("max_dl %.3g boot %.3g rel_beta %.3g rel_se %.3g df %s na/inf mismatches %d%s", dl, db, rel(5), rel(6), df_same, pat,
                      if (tag[[nm]] == "report") " (report-only case: rank-deficient / near-perfect, list any difference)" else ""))
})
write.table(do.call(rbind, rows), out, sep = "\t", quote = FALSE, row.names = FALSE)
