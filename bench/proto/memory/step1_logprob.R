# log density of the step-1 Stan model of a build at the baseline optimum (same function check)
# Rscript bench/proto/memory/step1_logprob.R <lib>
lib = commandArgs(TRUE)[1]
suppressPackageStartupMessages({ library(suez, lib.loc = lib); library(rstan) })
d = readRDS("bench/data/S.rds"); input = suez:::quantile_normalize(d$input); anno = d$anno
base = readRDS("bench/results/out/baseline_S_step1_c1_r1.rds")$optimized_fit$par
X = model.matrix(~ as.factor(anno$condition)); q = qr(X)
Q = qr.Q(q, complete = TRUE)[, -seq_len(q$rank), drop = FALSE]; resi = input %*% Q
same_ind = outer(anno$individual, anno$individual, "==") * 1
dat = list(YtY = crossprod(resi), N = ncol(resi), G = nrow(resi), P = 2L,
  x = list(diag(ncol(resi)), crossprod(Q, same_ind %*% Q)), P_adjustable = 10L)
m = suez:::stanmodels$suez_step_1
fit = sampling(m, data = dat, init = list(list(s = unname(base$s), s_adjustable = base$s_adjustable,
  x_adjustable = base$x_adjustable)), chains = 1, iter = 1, algorithm = "Fixed_param", refresh = 0)
pars = list(s = unname(base$s), s_adjustable = base$s_adjustable, x_adjustable = base$x_adjustable)
u = unconstrain_pars(fit, pars)
cat(sprintf("%s log_prob %.10f  |grad| %.3e\n", lib, log_prob(fit, u, adjust_transform = FALSE),
  sqrt(sum(grad_log_prob(fit, u, adjust_transform = FALSE)^2))))
