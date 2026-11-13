# Local-ancestry edge cases run end-to-end with one installed build: Rscript la_edge.R <lib dir> <out.rds>
# A: 3 ancestries, condition within individual, boot; B: 2 ancestries, ancestry 1 absent at some SNPs and
# rare (< 5 ALT) at others; C: q = 0 (single condition); D: 3 ancestries with ancestry 2 absent at every SNP of a
# gene and a covariate; E: between-individual condition (1 sample per individual), 2 ancestries.
a = commandArgs(TRUE); suppressPackageStartupMessages({ library(suez, lib.loc = a[1]); library(doMC) }); registerDoMC(1)
sim_la = function(ni, k, ns, theta) {
  p = matrix(runif(ns * k, 0.2, 0.6), ns, k); G = H = replicate(k, matrix(0, ns, ni), simplify = FALSE)
  for (h in 1:2) { an = sapply(1:ni, function(i) sample.int(k, 1, prob = theta[i, ]))
    for (j in 1:ns) { al = rbinom(ni, 1, p[cbind(j, an)]); for (b in 1:k) { H[[b]][j, ] = H[[b]][j, ] + (an == b); G[[b]][j, ] = G[[b]][j, ] + al * (an == b) } } }
  ids = paste0("I", 1:ni); sn = paste0("rs", 1:ns)
  f = function(L) { L = lapply(L, function(m) { dimnames(m) = list(sn, ids); m }); names(L) = paste0("anc", 1:k - 1); L }
  list(dosage = f(G), hapcount = f(H))
}
run = function(la, ni, reps, cond, ng = 3, covs = NULL) {
  ns = nrow(la$dosage[[1]]); ind = rep(paste0("I", 1:ni), each = reps); N = length(ind)
  g = la$dosage[[1]][1, ind]
  Y = t(sapply(1:ng, function(j) 0.4 * g * (1 + (j == 1) * cond) + rnorm(ni)[match(ind, paste0("I", 1:ni))] + rnorm(N)))
  dimnames(Y) = list(paste0("g", 1:ng), paste0("s", 1:N))
  snploc = data.frame(snpid = rownames(la$dosage[[1]]), chr = "1", pos = rep(1:ng, length.out = ns) * 1e6 + seq_len(ns))
  geneloc = data.frame(geneid = rownames(Y), chr = "1", left = 1:ng * 1e6, right = 1:ng * 1e6 + 10)
  set.seed(99); as.data.frame(map_interaction_qtl(Y, NULL, geneloc, snploc, data.frame(individual = ind, condition = cond),
    normalization_approach = "qq", permutation_approach = "boot", cisdist = 1e5, local_ancestry = la, covariates = covs))
}
set.seed(1); out = list()
ni = 150; th = matrix(rgamma(ni * 3, 2), ni); th = th / rowSums(th)
out$A_3anc = run(sim_la(ni, 3, 12, th), ni, 2, rep(0:1, ni))
la = sim_la(ni, 2, 12, cbind(rbeta(ni, 4, 2), 0)[, 1] |> (\(x) cbind(x, 1 - x))())
for (j in 1:3) { la$hapcount$anc0[j, ] = la$hapcount$anc0[j, ] + la$hapcount$anc1[j, ]; la$dosage$anc0[j, ] = la$dosage$anc0[j, ] + la$dosage$anc1[j, ]; la$hapcount$anc1[j, ] = 0; la$dosage$anc1[j, ] = 0 }
for (j in 4:6) { k = which(la$dosage$anc1[j, ] > 0); la$dosage$anc1[j, k[-(1:3)]] = 0 }   # 3 ALT carriers left on ancestry 1
out$B_absent_rare = run(la, ni, 2, rep(0:1, ni))
out$C_q0 = run(sim_la(ni, 2, 12, cbind(rep(0.6, ni), 0.4)), ni, 2, rep(1, 2 * ni))
la = sim_la(ni, 3, 12, th); la$hapcount$anc0 = la$hapcount$anc0 + la$hapcount$anc2; la$dosage$anc0 = la$dosage$anc0 + la$dosage$anc2
la$hapcount$anc2[] = 0; la$dosage$anc2[] = 0
out$D_3anc_one_absent = run(la, ni, 2, rep(0:1, ni), covs = cbind(age = rnorm(2 * ni)))
out$E_between = run(sim_la(200, 2, 12, cbind(rep(0.5, 200), 0.5)), 200, 1, rep(0:1, each = 100))
saveRDS(out, a[2])
