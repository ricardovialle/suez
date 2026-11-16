# Usage: Rscript bench/age_varying/sim_age.R bench/age_varying/results/sim_age.json, then python3 doc/src/build.py
# Linear (previous default) vs spline (varying-coefficient) age x genotype eQTL mapping in suez
suppressPackageStartupMessages({ library(suez); library(doMC); library(jsonlite) })
registerDoMC(8); set.seed(2026)

ni = 300                                              # individuals, 2 visits each, 4-10 years apart
age0 = runif(ni, 20, 85); age = c(age0, age0 + runif(ni, 4, 10))
ids = sprintf("I%03d", 1:ni); ind = rep(ids, 2); N = 2 * ni
a = (age - 55) / 20                                    # scaled age
shapes = list(                                         # beta(age) per ALT allele, mean ~ 0.35
  constant  = function(x) rep(0.35, length(x)),
  linear    = function(x) 0.35 + 0.25 * (x - 55) / 20,
  late      = function(x) 0.15 + 0.6 * pmax(0, (x - 60) / 20)^1.5,     # switches on after 60
  midlife   = function(x) 0.05 + 0.6 * exp(-((x - 55) / 12)^2)) # peaks at 55, reverses at the extremes
mu = function(x) 0.8 * sin((x - 20) / 25) - 0.3 * ((x - 55) / 20)^2   # nonlinear trajectory, independent of G

reps = 400; n_snp = 4
sim = function(shape) {
  G = matrix(rbinom(reps * n_snp * ni, 2, 0.3), reps * n_snp, ni, dimnames = list(paste0("rs", 1:(reps * n_snp)), ids))
  ui = rnorm(ni, 0, 0.5)
  Y = t(sapply(1:reps, function(g) {
    gs = G[(g - 1) * n_snp + 1, ind]
    mu(age) + shapes[[shape]](age) * gs + ui[match(ind, ids)] + rnorm(N, 0, 0.8)
  }))
  dimnames(Y) = list(paste0("g", 1:reps), paste0("s", 1:N))
  list(G = G, Y = Y,
       snploc = data.frame(snpid = rownames(G), chr = "1", pos = rep(1:reps, each = n_snp) * 1e6 + 1:n_snp),
       geneloc = data.frame(geneid = rownames(Y), chr = "1", left = 1:reps * 1e6, right = 1:reps * 1e6 + 10))
}
B = splines::ns(age, df = 3)
grid = seq(25, 90, 5); Bg = predict(B, grid); rownames(Bg) = grid
ag = cbind((grid - 55) / 20); rownames(ag) = grid
annoL = data.frame(individual = ind, condition = a)
annoS = data.frame(individual = ind); annoS$condition = B
run = function(d, anno, at) map_interaction_qtl(d$Y, d$G, d$geneloc, d$snploc, anno, normalization_approach = "none",
                                                 permutation_approach = "none", cisdist = 100, curve_at = at)

out = list(); t_l = t_s = 0
for (sh in names(shapes)) {
  d = sim(sh)
  t_l = t_l + system.time(rl <- run(d, annoL, ag))[3]; t_s = t_s + system.time(rs <- run(d, annoS, Bg))[3]
  causal = rownames(d$G)[(0:(reps - 1)) * n_snp + 1]
  cl = rl[rl$cis_snp %in% causal, ]; cs = rs[rs$cis_snp %in% causal, ]
  nl = rl[!rl$cis_snp %in% causal, ]; ns_ = rs[!rs$cis_snp %in% causal, ]   # non-causal SNPs: null for G x age
  # curves: median over replicates of estimates, and the first replicate with its band
  cv = function(r) { b = as.matrix(r[, paste0("beta_curve_", grid)]); s = as.matrix(r[, paste0("se_beta_curve_", grid)])
    list(median = unname(apply(b, 2, median)), ex = unname(b[1, ]), ex_se = unname(s[1, ]),
         rmse = sqrt(mean(sweep(b, 2, shapes[[sh]](grid))^2))) }
  out[[sh]] = list(
    truth = shapes[[sh]](grid),
    power_linear = mean(cl$p_interact < 0.05), power_spline = mean(cs$p_interact < 0.05),
    power5_linear = mean(cl$p_interact < 5e-5), power5_spline = mean(cs$p_interact < 5e-5),
    null_linear = mean(nl$p_interact < 0.05), null_spline = mean(ns_$p_interact < 0.05),
    lambda_linear = median(qchisq(nl$p_interact, 1, lower.tail = FALSE)) / qchisq(0.5, 1),
    lambda_spline = median(qchisq(ns_$p_interact, 1, lower.tail = FALSE)) / qchisq(0.5, 1),
    qq_linear = sort(-log10(rl$p_interact)), qq_spline = sort(-log10(rs$p_interact)),   # used for the constant shape (all SNPs null)
    curve_linear = cv(cl), curve_spline = cv(cs),
    mu_spline = unname(unlist(rs[1, paste0("mu_curve_", grid)])),
    mu_spline_se = unname(unlist(rs[1, paste0("se_mu_curve_", grid)])),
    mu_linear = unname(unlist(rl[1, paste0("mu_curve_", grid)])))
  cat(sh, "power L/S:", out[[sh]]$power_linear, out[[sh]]$power_spline, " null L/S:", out[[sh]]$null_linear, out[[sh]]$null_spline, "\n")
}
# true mu relative to the reference point (ns basis 0 at the lower boundary knot = min age; linear: age 55)
pc = list(linear = 10^-out$constant$qq_linear, spline = 10^-out$constant$qq_spline)   # constant shape: all SNPs null
calib = list(lam_linear = median(qchisq(pc$linear, 1, lower.tail = FALSE)) / qchisq(0.5, 1),
             lam_spline = median(qchisq(pc$spline, 1, lower.tail = FALSE)) / qchisq(0.5, 1),
             fp_linear = mean(pc$linear < 0.05), fp_spline = mean(pc$spline < 0.05))
res = list(grid = grid, shapes = out, calib = calib, mu_true = mu(grid), age_min = min(age), n_ind = ni, n_samples = N, reps = reps,
           n_snp = n_snp, sec_linear = unname(t_l), sec_spline = unname(t_s),
           age_hist = hist(age, breaks = seq(15, 100, 5), plot = FALSE)$counts)
writeLines(toJSON(res, digits = 5, auto_unbox = TRUE), commandArgs(TRUE)[1])
