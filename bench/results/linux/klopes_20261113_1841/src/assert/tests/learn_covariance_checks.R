# Checks of learn_covariance (suez step 1): variance components are recovered from data
# simulated with known noise, same-individual and latent-factor variance, for a factor and a
# continuous condition with large condition effects per gene.
suppressPackageStartupMessages({ library(suez); library(doMC) })
fails = 0
check = function(name, ok, detail = "") {
  cat(if (isTRUE(ok)) "PASS" else "FAIL", name, detail, "\n")
  if (!isTRUE(ok)) fails <<- fails + 1
}
simulate = function(condition, ind, G = 400, s_noise = 0.5, s_ind = 0.3, s_f = 2) {
  N = length(ind)
  f = rnorm(N); f = f / sqrt(sum(f^2))
  Sigma = s_noise * diag(N) + s_ind * outer(ind, ind, "==") + s_f * tcrossprod(f)
  X = model.matrix(~condition)
  E = matrix(rnorm(G * N), G) %*% chol(Sigma)
  B = matrix(rnorm(G * ncol(X), sd = 2), G)
  list(y = B %*% t(X) + E, f = f)
}
set.seed(3)
ind = rep(sprintf("I%02d", 1:40), each = 2)
for (case in c("factor", "continuous")) {
  cond = if (case == "factor") factor(rep(c("ctl", "stim"), 40)) else rnorm(80)
  d = simulate(cond, ind)
  fit = learn_covariance(d$y, data.frame(individual = ind, condition = cond), P = 2)
  s = fit$optimized_fit$par$s
  sf = max(fit$optimized_fit$par$s_adjustable)
  check(paste(case, "condition: converged"), fit$optimized_fit$return_code == 0)
  check(paste(case, "condition: noise and same-individual variances"),
        abs(s["isotropic_noise"] / 0.5 - 1) < 0.1 && abs(s["same_ind"] / 0.3 - 1) < 0.25,
        sprintf("(noise %.3f vs 0.5, same_ind %.3f vs 0.3, factor %.2f vs 2)", s["isotropic_noise"], s["same_ind"], sf))
  check(paste(case, "condition: cov_suez is a PSD samples x samples matrix"),
        all(dim(fit$cov_suez) == 80) && min(eigen(fit$cov_suez, TRUE, TRUE)$values) > -1e-8)
}
cat(fails, "failure(s)\n")
if (fails > 0) stop(fails, " check(s) failed")
