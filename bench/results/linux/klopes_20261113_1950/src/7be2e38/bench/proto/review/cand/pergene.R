# per-gene l0 error vs exact and the fitted x = log10(s2k/s2) for the extreme cases (p7 vs p3)
src = readLines("bench/proto/review/cand/extreme_cases.R"); eval(parse(text = src[1:grep("^options", src)]))
for (nm in c("default_h_to_1", "grm_reps_x1e6_h_to_1", "null_noise_1em13")) { P = cases[[nm]]
  for (k in names(impl)) { r = run_impl(impl[[k]], P); l0 = sapply(r$out, `[[`, "l0")
    ex = sapply(P$genes, function(g) exact(P$UX0, g$Uy, P$lambda))
    cat(sprintf("%-22s %-17s l0 - exact per gene: %s\n", nm, k, paste(sprintf("%.2g", l0 - ex), collapse = " "))) } }
# where is the optimum? x of the exact fit for default_h_to_1
P = cases$default_h_to_1; lam = P$lambda
for (g in P$genes) { f = function(x) { h = 1/(1+10^-x); d = h*lam+1-h; W = 1/d; b = qr.coef(qr(P$UX0*sqrt(W)), g$Uy*sqrt(W)); r = g$Uy - P$UX0 %*% b; -0.5*sum(log(sum(r^2*W)/length(r)*d)) }
  cat("argmax x:", optimize(f, c(-12, 12), maximum = TRUE)$maximum, "\n") }
