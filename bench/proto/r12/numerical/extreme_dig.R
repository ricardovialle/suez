# Dig into extreme-case differences vs base12: per-SNP exact l_geno, and base12's own floor (Uy x (1 + 1e-15)).
#   ONLY=<case,...> Rscript bench/proto/r12/numerical/extreme_dig.R <proto.cpp> ...
Sys.setenv(SUEZ_HARNESS_LIB = "bench/lib/base12")
src = readLines("bench/proto/review/cand/extreme_cases.R")
eval(parse(text = src[grep("^source", src):(grep("^options", src) - 1)]))   # harness, impl list, mk, exact, cases
for (nm in names(cases)) {
  P = cases[[nm]]; b = run_impl(baseline_fit, P)
  Pp = P; Pp$genes = lapply(P$genes, function(g) { g$Uy = g$Uy * (1 + 1e-15); g }); bp = run_impl(baseline_fit, Pp)
  fl = compare(b, bp)
  cat(sprintf("\n## %s: base12 floor (Uy x (1+1e-15)): dl_geno %.1e dl_int %.1e dl_boot %.1e/%.1e\n", nm, fl["l_geno"], fl["l_interact"], fl["l_boot_geno"], fl["l_boot_interact"]))
  for (k in names(impl)[-1]) {
    r = run_impl(impl[[k]], P); cm = compare(b, r)
    B = do.call(rbind, lapply(b$out, `[[`, "res")); R = do.call(rbind, lapply(r$out, `[[`, "res"))
    i = which.max(abs(B[, 1] - R[, 1])); gi = (i - 1) %/% ncol(P$genes[[1]]$G) + 1; si = (i - 1) %% ncol(P$genes[[1]]$G) + 1
    g = P$genes[[gi]]; ex = exact(cbind(P$UX0, crossprod(P$U, std(g$G[, si]))), g$Uy, P$lambda)
    cat(sprintf("  %s: dl_geno %.1e dl_int %.1e dl_boot %.1e/%.1e; worst SNP (gene %d, snp %d): base - exact %.2e, proto - exact %.2e\n", k,
                cm["l_geno"], cm["l_interact"], cm["l_boot_geno"], cm["l_boot_interact"], gi, si, B[i, 1] - ex, R[i, 1] - ex))
  }
}
