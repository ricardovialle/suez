# Checks of map_interaction_qtl(local_ancestry = ...) (run by R CMD check, or Rscript tests/local_ancestry_checks.R):
#   1. every likelihood equals the exact profile-ML reference (2 conditions, repeated samples, covariate)
#   2. one sample per individual and the default kernel: ancestry betas, SEs and p-values equal Tractor's
#      lm(y ~ cov + LA + G_1..G_k)
#   3. three ancestries; an ancestry with < 5 ALT alleles is dropped; read_tractor() reads Tractor's files
suppressPackageStartupMessages({ library(suez); library(doMC) })
registerDoMC(1)
fails = 0
check = function(name, ok, detail = "") {
  cat(if (isTRUE(ok)) "PASS" else "FAIL", name, detail, "\n")
  if (!isTRUE(ok)) fails <<- fails + 1
}

sim_la = function(ni, k, ns, theta) {   # haplotype ancestry and alleles, no LD
  p = matrix(runif(ns * k, 0.2, 0.6), ns, k)
  G = H = replicate(k, matrix(0, ns, ni), simplify = FALSE)
  for (h in 1:2) {
    a = sapply(1:ni, function(i) sample.int(k, 1, prob = theta[i, ]))
    for (j in 1:ns) {
      al = rbinom(ni, 1, p[cbind(j, a)])
      for (b in 1:k) { H[[b]][j, ] = H[[b]][j, ] + (a == b); G[[b]][j, ] = G[[b]][j, ] + al * (a == b) }
    }
  }
  nm = paste0("anc", 1:k - 1)
  ids = paste0("I", 1:ni); sn = paste0("rs", 1:ns)
  f = function(L) { L = lapply(L, function(m) { dimnames(m) = list(sn, ids); m }); names(L) = nm; L }
  list(dosage = f(G), hapcount = f(H))
}
exact_ml = function(UX, Uy, lam) {
  UX = UX[, colSums(UX^2) > 1e-12, drop = FALSE]
  f = function(x) { h = 1 / (1 + 10^-x); d = h * lam + (1 - h); W = 1 / d
    b = solve(crossprod(UX, UX * W), crossprod(UX, Uy * W))
    -0.5 * length(Uy) - 0.5 * sum(log(sum((Uy - UX %*% b)^2 * W) / length(Uy) * d)) }
  xs = seq(-8, 8, length.out = 801); v = sapply(xs, f); pk = which(diff(sign(diff(v))) == -2) + 1
  max(v, sapply(pk, function(i) optimize(f, xs[c(i - 1, i + 1)], maximum = TRUE, tol = 1e-10)$objective))
}
std = function(x) if (sd(x) > 0) (x - mean(x)) / sd(x) else 0 * x

# 1. two conditions per individual, 2 ancestries, a covariate
set.seed(21)
ni = 150; k = 2; ns = 6
theta = cbind(rbeta(ni, 4, 2)); theta = cbind(theta, 1 - theta)
la = sim_la(ni, k, ns, theta)
ind = rep(paste0("I", 1:ni), 2); cond = rep(0:1, each = ni); N = 2 * ni
g0 = la$dosage$anc0[1, ind]; g1 = la$dosage$anc1[1, ind]; lh = la$hapcount$anc0[1, ind]
y = 0.3 * cond + 0.5 * g0 + 0.8 * g0 * cond + 0.2 * lh + theta[match(ind, paste0("I", 1:ni)), 1] +
  rnorm(ni, 0, 0.6)[match(ind, paste0("I", 1:ni))] + rnorm(N, 0, 0.7)
Y = rbind(gA = y); colnames(Y) = paste0("s", 1:N)
snploc = data.frame(snpid = rownames(la$dosage$anc0), chr = "1", pos = 1000 + 1:ns)
geneloc = data.frame(geneid = "gA", chr = "1", left = 1000, right = 1001)
anno = data.frame(individual = ind, condition = cond)
cov1 = cbind(global = theta[match(ind, paste0("I", 1:ni)), 1])
set.seed(1)
r = map_interaction_qtl(Y, NULL, geneloc, snploc, anno, normalization_approach = "none", cisdist = 100,
                        local_ancestry = la, covariates = cov1, permutation_approach = "boot")
check("one row per SNP, all tests", nrow(r) == ns && !anyNA(r[, c("p_geno", "p_het", "p_interact", "p_interact_het")]))
K = outer(ind, ind, "==") * 1
e = suez:::kernel_contrast_eigen(K); U = e$vectors; lam = e$values
UX0 = crossprod(U, cbind(std(cond), std(cov1[, 1])))
worst = 0
for (i in seq_len(ns)) {
  s = r$cis_snp[i]; l_ = la$hapcount$anc0[s, ind]; a0 = la$dosage$anc0[s, ind]; a1 = la$dosage$anc1[s, ind]
  cols = cbind(l_, a0, a1, a0 + a1, l_ * cond, a0 * cond, a1 * cond, (a0 + a1) * cond)
  cols = crossprod(U, apply(cols, 2, std))
  mods = list(l0 = 1, l_geno = 1:3, l_geno_shared = c(1, 4), l_interact_null = c(1:3, 5),
              l_interact = c(1:3, 5:7), l_interact_shared = c(1:3, 5, 8))
  for (m in names(mods)) worst = max(worst, abs(r[i, m] - exact_ml(cbind(UX0, cols[, mods[[m]], drop = FALSE]), drop(crossprod(U, y)), lam)))
}
check("local-ancestry log-likelihoods = exact profile ML", worst < 1e-6, sprintf("(%d fits, max diff %.2g)", 6 * ns, worst))
check("bootstrap columns", all(c("l_boot_interact_null", "l_boot_interact") %in% names(r)) && !anyNA(r$l_boot_interact))
check("G0 x condition detected", r$p_interact[1] < 1e-3 && r$p_interact_anc0[1] < r$p_interact_anc1[1])

# 2. one sample per individual, single condition: Tractor's linear model
set.seed(22)
ni = 120; la2 = sim_la(ni, 2, 5, cbind(rep(0.7, ni), rep(0.3, ni)))
ids = paste0("I", 1:ni); age = rnorm(ni)
y2 = 0.4 * la2$dosage$anc0[2, ] + 0.3 * age + rnorm(ni)
Y2 = rbind(gB = y2); colnames(Y2) = ids
r2 = map_interaction_qtl(Y2, NULL, data.frame(geneid = "gB", chr = "1", left = 1000, right = 1001), snploc[1:5, ],
                         data.frame(individual = ids, condition = 1), normalization_approach = "none", cisdist = 100,
                         local_ancestry = la2, covariates = cbind(age = age), permutation_approach = "none")
db = ds = dp = 0
for (i in seq_len(nrow(r2))) {
  s = r2$cis_snp[i]
  fit = summary(lm(y2 ~ age + la2$hapcount$anc0[s, ] + la2$dosage$anc0[s, ] + la2$dosage$anc1[s, ]))$coefficients
  db = max(db, abs(c(r2$beta_anc0[i], r2$beta_anc1[i]) - fit[4:5, 1]))
  ds = max(ds, abs(c(r2$se_anc0[i], r2$se_anc1[i]) - fit[4:5, 2]))
  dp = max(dp, abs(log(c(r2$p_anc0[i], r2$p_anc1[i]) / fit[4:5, 4])))
}
check("one sample per individual: betas = Tractor lm", db < 1e-8, sprintf("(max diff %.2g)", db))
check("one sample per individual: SEs and p-values = Tractor lm", ds < 1e-8 && dp < 1e-8, sprintf("(max diff %.2g, %.2g)", ds, dp))
check("no interaction columns for a single condition", all(is.na(r2$p_interact)) && all(r2$df_interact == 0))

# 3. three ancestries, a rare ancestry-specific allele, read_tractor()
set.seed(23)
ni = 200; th3 = matrix(rep(c(0.5, 0.3, 0.2), each = ni), ni)
la3 = sim_la(ni, 3, 4, th3)
la3$dosage$anc2[4, ] = 0; la3$dosage$anc2[4, which(la3$hapcount$anc2[4, ] > 0)[1:2]] = 1   # 2 ALT alleles on anc2
la3$dosage$anc0[3, ] = la3$hapcount$anc0[3, ]                                                # every anc0 haplotype ALT
y3 = 1.0 * la3$dosage$anc1[1, ] + rnorm(ni); Y3 = rbind(gC = y3); colnames(Y3) = paste0("I", 1:ni)
r3 = map_interaction_qtl(Y3, NULL, data.frame(geneid = "gC", chr = "1", left = 1000, right = 1001), snploc[1:4, ],
                         data.frame(individual = paste0("I", 1:ni), condition = 1), normalization_approach = "none",
                         cisdist = 100, local_ancestry = la3, permutation_approach = "none")
check("three ancestries: joint df 3, heterogeneity df 2", r3$df_geno[1] == 3 && r3$df_het[1] == 2 && r3$p_geno[1] < 1e-3)
check("ancestry with < 5 ALT alleles dropped", r3$df_geno[4] == 2 && is.na(r3$beta_anc2[4]) && !is.na(r3$beta_anc0[4]))
check("ancestry with < 5 REF alleles dropped (collinear with LA)", r3$df_geno[3] == 2 && is.na(r3$beta_anc0[3]) && is.finite(r3$p_geno[3]) && r3$p_geno[3] > 0)

# third ancestry absent at the locus: the local-ancestry columns are collinear (hapcounts sum to 2); the
# genotype effects stay identified and equal those of the two-ancestry model
la4 = sim_la(ni, 2, 4, cbind(rep(0.6, ni), rep(0.4, ni)))
for (m in c("dosage", "hapcount")) { la4[[m]]$anc2 = la4[[m]]$anc0 * 0 }
y4 = 0.8 * la4$dosage$anc0[1, ] + rnorm(ni); Y4 = rbind(gD = y4); colnames(Y4) = paste0("I", 1:ni)
fit4 = function(la) map_interaction_qtl(Y4, NULL, data.frame(geneid = "gD", chr = "1", left = 1000, right = 1001), snploc[1:4, ],
  data.frame(individual = paste0("I", 1:ni), condition = 1), normalization_approach = "none", cisdist = 100,
  local_ancestry = la, permutation_approach = "none")
r4 = fit4(la4); r2a = fit4(list(dosage = la4$dosage[1:2], hapcount = la4$hapcount[1:2]))
check("absent third ancestry: identified genotype effects equal the two-ancestry fit",
      all(is.finite(r4$se_anc0)) && all(r4$se_anc0 > 0) && max(abs(r4$beta_anc0 - r2a$beta_anc0)) < 1e-8 &&
        max(abs(r4$l_geno - r2a$l_geno)) < 1e-8, sprintf("(max |dbeta| %.2g)", max(abs(r4$beta_anc0 - r2a$beta_anc0))))

td = file.path(tempdir(), "tractor_files"); dir.create(td, showWarnings = FALSE)
for (a in 0:1) for (type in c("dosage", "hapcount")) {
  m = la2[[type]][[a + 1]]
  out = data.frame(CHROM = "chr22", POS = snploc$pos[1:5], ID = c(".", "rs2", "rs3", ".", "rs5"), REF = "A", ALT = "G", m, check.names = FALSE)
  write.table(out, gzfile(file.path(td, sprintf("toy.anc%d.%s.txt.gz", a, type))), sep = "\t", quote = FALSE, row.names = FALSE)
}
rt = read_tractor(file.path(td, "toy"), ancestry_names = c("AFR", "EUR"))
check("read_tractor reads Tractor files", identical(names(rt$dosage), c("AFR", "EUR")) &&
        all(rt$dosage$EUR == la2$dosage$anc1) && all(rt$hapcount$AFR == la2$hapcount$anc0) &&
        identical(rt$snploc$snpid, c("chr22:1001:A:G", "rs2", "rs3", "chr22:1004:A:G", "rs5")))
unlink(td, recursive = TRUE)

cat(fails, "failure(s)\n")
if (fails > 0) stop(fails, " check(s) failed")
