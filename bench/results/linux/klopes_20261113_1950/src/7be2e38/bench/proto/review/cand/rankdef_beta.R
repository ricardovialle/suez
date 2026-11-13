# beta/se of an identified term when X0 is rank-deficient (duplicated covariate): candidate vs baseline vs the
# same fit with the duplicate removed (the identified answer).   Rscript rankdef_beta.R cand.cpp
source("bench/proto/harness.R"); e = new.env(); Rcpp::sourceCpp(commandArgs(TRUE)[1], env = e)
set.seed(4); ni = 100; ind = rep(1:ni, each = 2); N = 2 * ni; cnd = rep(c(0, 1), ni)
K = outer(ind, ind, "==") * 1; eig = suez:::kernel_contrast_eigen(K); U = eig$vectors
std = function(x) (x - mean(x)) / sd(x); cv = rnorm(N)
G = sapply(runif(10, 0.1, 0.5), function(p) rbinom(ni, 2, p))[ind, ]; y = 0.3 * cnd + 0.3 * G[, 1] + rnorm(ni)[ind] + rnorm(N)
UXd = crossprod(U, cbind(std(cnd), std(cv), std(cv))); UXr = crossprod(U, cbind(std(cnd), std(cv))); Uy = drop(crossprod(U, y))
ref = baseline_fit(U, eig$values, Uy, UXr, G, cbind(cnd), FALSE)$res
for (nm in c("baseline", "candidate")) { f = if (nm == "baseline") baseline_fit else e$fit_gene_proto
  r = f(U, eig$values, Uy, UXd, G, cbind(cnd), FALSE)$res
  cat(sprintf("%-9s max|dl| %.2g  beta: max rel err %.2g, NaN %d  se: max rel err %.2g, NaN %d, zero %d\n", nm, max(abs(r[, 1:2] - ref[, 1:2])),
    max(abs(r[, 4] - ref[, 4]) / ref[, 5], na.rm = TRUE), sum(is.na(r[, 4])), max(abs(r[, 5] / ref[, 5] - 1), na.rm = TRUE), sum(is.na(r[, 5])), sum(r[, 5] == 0, na.rm = TRUE))) }
