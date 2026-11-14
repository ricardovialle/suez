# Which entries differ in NA pattern on edge case collinear_interaction_and_snps (p_prof2 vs base12)
Sys.setenv(SUEZ_HARNESS_LIB = "bench/lib/base12")
src = readLines("bench/proto/review/edge_cases.R"); src = src[seq_len(grep("^options\\(width = 200\\)", src)[1] - 1)]
src = sub("a = commandArgs(TRUE)", "a = 'lib:'", src, fixed = TRUE); eval(parse(text = src))
e = new.env(); Rcpp::sourceCpp("bench/proto/r12/numerical/p_prof2.cpp", env = e)
P = cases$collinear_interaction_and_snps; b = run_impl(baseline_fit, P); r = run_impl(e$fit_gene_proto, P)
B = do.call(rbind, lapply(b$out, `[[`, "res")); R = do.call(rbind, lapply(r$out, `[[`, "res"))
w = which(is.na(B) != is.na(R), arr.ind = TRUE); print(w); for (i in unique(w[, 1])) { print(rbind(base = B[i, ], cand = R[i, ])) }
# base12's own sensitivity: Uy x (1 + k 1e-15)
for (k in 1:5) { Pp = P; Pp$genes = lapply(P$genes, function(g) { g$Uy = g$Uy * (1 + k * 1e-15); g }); bp = run_impl(baseline_fit, Pp)
  Bp = do.call(rbind, lapply(bp$out, `[[`, "res")); cat("perturb", k, "NA mismatches vs base:", sum(is.na(B) != is.na(Bp)), " SNP12 beta:", Bp[12, 4], "\n") }
