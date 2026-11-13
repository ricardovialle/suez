# Interleaved timing of kernel-level implementations (fairer than back-to-back run_impl() calls while other
# jobs share the machine): one warm-up, then reps rounds, each running every implementation once, in an
# order rotated every round. Reports per-round elapsed and CPU (user+sys of this process) ratios vs the first.
#   source("bench/proto/harness.R"); source("bench/proto/review/interleave.R")
#   Rcpp::sourceCpp("bench/proto/p0_baseline.cpp"); p0 = fit_gene_proto
#   Rcpp::sourceCpp("cand.cpp"); time_interleaved(list(p0 = p0, cand = fit_gene_proto), prepare("M_grm", 30), reps = 5)
# Note: sourceCpp'd prototypes must be timed against p0 compiled the same way, not against the installed .Call.
time_interleaved = function(impls, P, reps = 5) {
  one = function(f) { t0 = proc.time(); for (g in P$genes) { set.seed(g$seed); f(P$U, P$lambda, g$Uy, P$UX0, g$G, P$C, P$boot) }
                      t = proc.time() - t0; c(wall = t[["elapsed"]], cpu = t[["user.self"]] + t[["sys.self"]]) }
  for (f in impls) one(f)                                           # warm-up (page faults, caches, lazy loads)
  k = length(impls); W = C = matrix(NA, reps, k, dimnames = list(NULL, names(impls)))
  for (r in seq_len(reps)) for (j in (seq_len(k) + r - 2) %% k + 1) { x = one(impls[[j]]); W[r, j] = x[["wall"]]; C[r, j] = x[["cpu"]] }
  out = data.frame(impl = names(impls), median_wall = apply(W, 2, median), median_cpu = apply(C, 2, median),
                   ratio_wall = apply(W / W[, 1], 2, median), ratio_cpu = apply(C / C[, 1], 2, median),
                   ratio_wall_min = apply(W / W[, 1], 2, min), ratio_wall_max = apply(W / W[, 1], 2, max))
  print(out, row.names = FALSE, digits = 3); invisible(list(summary = out, wall = W, cpu = C))
}
