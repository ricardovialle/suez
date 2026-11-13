source("bench/proto/harness.R"); e = new.env(); Rcpp::sourceCpp("bench/proto/numerical/dbg/p9_dbg.cpp", env = e)
set.seed(5); ni = 100; ind = rep(sprintf("i%03d", 1:ni), each = 2); N = 2 * ni; cnd = rep(c("a", "b"), ni)
Kd = outer(ind, ind, "==") * 1; eig = suez:::kernel_contrast_eigen(Kd)
X0 = suez:::standardize_cols(model.matrix(~ c, data.frame(c = cnd))[, -1, drop = FALSE]); UX = crossprod(eig$vectors, X0)
u = rnorm(ni)[rep(1:ni, each = 2)]; y = 0.3 * (cnd == "b") + u + 1e-4 * rnorm(N); Uy = drop(crossprod(eig$vectors, y)); lam = eig$values
f = function(h) { d = h * lam + (1 - h); W = 1 / d; b = qr.coef(qr(UX * sqrt(W)), Uy * sqrt(W)); r = Uy - UX %*% b
  c(-0.5 * length(Uy) - 0.5 * sum(log(sum(r^2 * W) / length(Uy) * d)), sum(log(d)), sum((lam - 1) * W), sum(((lam - 1) * W)^2)) }
for (x in c(0, 3, 6, 8, 10, 11, 12)) { h = 1 / (1 + 10^-x); a = e$prof_dbg(lam, Uy, UX, h); ex = f(h)
  cat(sprintf("x=%4.1f null_ok=%d R=%d n0=%d  dv=%.3g dLD=%.3g dT1=%.3g dT2=%.3g\n", x, a[6], a[4], a[5], a[1] - ex[1], a[7] - ex[2], a[8] - ex[3], a[9] - ex[4])) }
