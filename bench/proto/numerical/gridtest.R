# direct comparison of the grid values (data V1, V2 and bootstrap W0, W1, W2) of p5 vs p0, over all genes
setwd("/Users/rvialle/suez"); source("bench/proto/harness.R")
ld = function(f) { e = new.env(); Rcpp::sourceCpp(f, env = e, cacheDir = "bench/proto/numerical/.cache"); e }
a = ld("bench/proto/numerical/gridtest_p0.cpp"); b = ld("bench/proto/numerical/gridtest_p5.cpp")
for (d in commandArgs(TRUE)) { P = prepare(d, 20); worst = 0; peakdiff = 0; ninf = 0
  for (g in P$genes) { set.seed(g$seed); a$fit_gene_proto(P$U, P$lambda, g$Uy, P$UX0, g$G, P$C, P$boot); A = a$last_grids()
    set.seed(g$seed); b$fit_gene_proto(P$U, P$lambda, g$Uy, P$UX0, g$G, P$C, P$boot); B = b$last_grids()
    for (j in seq_len(if (ncol(P$C)) 2 else 1)) { x = A[[j]]; y = B[[j]]; ok = is.finite(x) & is.finite(y); ninf = ninf + sum(is.finite(x) != is.finite(y))
      dd = ifelse(ok, abs(x - y), 0); if (max(dd) > worst) { worst = max(dd); wi = which(dd == worst, arr.ind = TRUE)[1, ]; wv = x[wi[1], wi[2]]; wj = j }; peakdiff = peakdiff + sum(apply(x, 1, which.max) != apply(y, 1, which.max)) } }
  cat(sprintf("%s: grid max|diff|=%.2e  argmax mismatches=%d  finiteness mismatches=%d  fallbacks=%g  worst at grid col %d (h=%.6g) model %d value %.4g\n", d, worst, peakdiff, ninf, b$fallbacks(), wi[2], c(0, 1/(1+10^-seq(-5,5,length.out=80)), 1)[wi[2]], wj, wv)) }
