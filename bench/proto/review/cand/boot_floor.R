# Bootstrap noise floor vs n: baseline vs baseline on Uy * (1 + 1e-15 N(0,1)), and candidates vs baseline;
# l_boot columns and the bootstrap interaction LR 2 (l_boot_interact - l_boot_geno).
source("bench/proto/harness.R"); files = commandArgs(TRUE); impl = list()
for (f in files) { e = new.env(); Rcpp::sourceCpp(f, env = e); impl[[sub("\\.cpp$", "", basename(f))]] = e$fit_gene_proto }
res = function(r) do.call(rbind, lapply(r$out, `[[`, "res"))
for (ds in c("S", "M_grm", "L_grm", "XL")) {
  P = prepare(ds, c(S = 30, M_grm = 20, L_grm = 8, XL = 3)[[ds]]); A = res(run_impl(baseline_fit, P))
  Pf = P; set.seed(3); for (i in seq_along(Pf$genes)) Pf$genes[[i]]$Uy = Pf$genes[[i]]$Uy * (1 + 1e-15 * rnorm(length(Pf$genes[[i]]$Uy)))
  B = list(floor = res(run_impl(baseline_fit, Pf))); for (k in names(impl)) B[[k]] = res(run_impl(impl[[k]], P))
  for (k in names(B)) { X = B[[k]]
    cat(sprintf("%-6s n=%4d %-17s max|dl_boot| %.2g  max|d boot LR| %.2g  max|d data LR| %.2g\n", ds, length(P$lambda), k,
      max(abs(A[, 6:7] - X[, 6:7]), na.rm = TRUE), max(abs(2 * (A[, 7] - A[, 6]) - 2 * (X[, 7] - X[, 6])), na.rm = TRUE),
      max(abs(2 * (A[, 2] - A[, 1]) - 2 * (X[, 2] - X[, 1])), na.rm = TRUE))) }
}
