# Inspection of the baseline on the edge cases: convergence vs an exact optimizer, kernel-scale invariance,
# and what it returns on rank-deficient / near-perfect designs.   Rscript bench/proto/review/baseline_checks.R
source("bench/proto/review/edge_cases.R")   # builds `cases`; also runs the suite on p0 (cheap)
ll_exact = function(UX, Uy, lam) max(exact_peaks(UX, Uy, lam, ends = TRUE))
for (nm in c("default_kernel_balanced", "grm_numeric_cond_h0_to_1", "wide_spectrum_kernel_x1", "bimodal_likelihood", "small_odd_n31")) {
  P = cases[[nm]]; b = run_impl(baseline_fit, P)
  d = unlist(lapply(seq_along(P$genes), function(i) { g = P$genes[[i]]
    e0 = if (ncol(P$UX0)) ll_exact(P$UX0, g$Uy, P$lambda) else NA
    e1 = sapply(1:min(5, ncol(g$G)), function(s) ll_exact(cbind(P$UX0, crossprod(P$U, std(g$G[, s]))), g$Uy, P$lambda))
    c(b$out[[i]]$l0 - e0, b$out[[i]]$res[1:length(e1), 1] - e1) }))
  cat(sprintf("%-28s baseline - exact optimum: min %.2g max %.2g (negative = baseline below the optimum)\n", nm, min(d, na.rm = TRUE), max(d, na.rm = TRUE)))
}
r = lapply(c("wide_spectrum_kernel_x1", "wide_spectrum_kernel_x1e+06", "wide_spectrum_kernel_x1e-06"), function(nm) run_impl(baseline_fit, cases[[nm]]))
cat("kernel scale invariance (x1e6, x1e-6 vs x1):\n"); print(rbind(compare(r[[1]], r[[2]]), compare(r[[1]], r[[3]])))
for (nm in c("collinear_interaction_and_snps", "duplicated_covariate", "near_perfect_fit")) {
  b = run_impl(baseline_fit, cases[[nm]]); cat("\n", nm, "(gene 1; l0 =", b$out[[1]]$l0, ")\n"); print(signif(b$out[[1]]$res[1:min(6, nrow(b$out[[1]]$res)), ], 6))
}
