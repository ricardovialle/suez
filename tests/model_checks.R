# Checks of map_interaction_qtl's model fits (run by R CMD check, or Rscript tests/model_checks.R):
#   1. every likelihood-ratio test equals the exact profile-ML LRT (grid + optimize over
#      h = s2k / (s2k + s2), b and the scale profiled out) on the intercept-free contrasts, with
#      repeated samples, a centred GRM (K 1 = 0) and a condition varying within individuals,
#      including genes whose ML sits on the s2 = 0 or s2k = 0 boundary
#   2. GLS effect sizes, bootstrap columns
#   3. errors: non-PSD kernel, unidentifiable residual variance
suppressPackageStartupMessages({ library(suez); library(doMC) })
registerDoMC(1)
fails = 0
check = function(name, ok, detail = "") {
  cat(if (isTRUE(ok)) "PASS" else "FAIL", name, detail, "\n")
  if (!isTRUE(ok)) fails <<- fails + 1
}
errmsg = function(expr) tryCatch({ expr; "" }, error = function(e) conditionMessage(e))

profile_ml = function(UX, Uy, lam) {
  f = function(h) {
    d = h * lam + (1 - h); if (any(d <= 0)) return(-Inf); W = 1 / d
    b = if (ncol(UX)) solve(crossprod(UX, UX * W), crossprod(UX, Uy * W)) else numeric(0)
    R = Uy - UX %*% b; s = sum(R^2 * W) / length(Uy)
    -0.5 * length(Uy) - 0.5 * sum(log(s * d))
  }
  hs = seq(0, 1, length.out = 201); v = sapply(hs, f); i = which.max(v)
  max(v[i], optimize(f, hs[c(max(1, i - 1), min(201, i + 1))], maximum = TRUE, tol = 1e-12)$objective)
}
std = function(x) (x - mean(x)) / sd(x)

set.seed(7)
nd = 120; reps = sample(1:2, nd, TRUE, prob = c(0.5, 0.5))
inds = sprintf("I%03d", 1:nd); ind = rep(inds, reps); N = length(ind)
A = matrix(rbinom(nd * 2000, 2, 0.3), nd)
Z = scale(A); grm = tcrossprod(Z) / ncol(Z); dimnames(grm) = list(inds, inds)
K = grm[ind, ind]
check("centred donor GRM is singular along 1", abs(sum(grm %*% rep(1, nd))) < 1e-6)
cond = rnorm(N)                                         # varies within individuals
n_gene = 8; n_snp = 6
G = matrix(rbinom(n_gene * n_snp * nd, 2, 0.3), n_gene * n_snp, nd,
           dimnames = list(paste0("rs", 1:(n_gene * n_snp)), inds))
snploc = data.frame(snpid = rownames(G), chr = "1", pos = rep(1:n_gene, each = n_snp) * 1e6 + rep(1:n_snp, n_gene))
geneloc = data.frame(geneid = paste0("g", 1:n_gene), chr = "1", left = 1:n_gene * 1e6, right = 1:n_gene * 1e6 + 10)
u = drop(t(chol(K + 1e-6 * diag(N))) %*% rnorm(N))
Y = t(sapply(1:n_gene, function(g) {
  h = c(0, 0.3, 0.8, 1)[g %% 4 + 1]                     # includes both boundaries
  gs = G[(g - 1) * n_snp + 1, ind]
  0.3 * cond + 0.4 * gs + 0.3 * gs * cond + sqrt(h) * std(u) + sqrt(1 - h) * rnorm(N)
}))
dimnames(Y) = list(geneloc$geneid, paste0("s", 1:N))
anno = data.frame(individual = ind, condition = cond)
res = map_interaction_qtl(Y, G, geneloc, snploc, anno, sample_kernel = K, normalization_approach = "none",
                          permutation_approach = "boot", cisdist = 100)

check("one row per gene x SNP, no failed fits",
      nrow(res) == n_gene * n_snp && !anyNA(res[, c("l_geno", "l_interact", "p_geno", "p_interact")]))
Q = qr.Q(qr(matrix(1, N, 1)), complete = TRUE)[, -1]
e = eigen(crossprod(Q, K %*% Q), symmetric = TRUE); U = Q %*% e$vectors; lam = pmax(e$values, 0)
Uy = function(g) drop(crossprod(U, Y[g, ]))
UX0 = crossprod(U, cbind(std(cond)))
worst = 0
for (i in seq_len(nrow(res))) {
  g = res$gene[i]; gs = G[res$cis_snp[i], ind]
  UX1 = cbind(UX0, crossprod(U, std(gs))); UX2 = cbind(UX1, crossprod(U, std(gs * cond)))
  p0 = profile_ml(UX0, Uy(g), lam); p1 = profile_ml(UX1, Uy(g), lam); p2 = profile_ml(UX2, Uy(g), lam)
  d = c(2 * (res$l_geno[i] - res$l0[i]) - 2 * (p1 - p0), 2 * (res$l_interact[i] - res$l_geno[i]) - 2 * (p2 - p1))
  worst = max(worst, abs(d))
}
check("LRT = exact profile-ML LRT", worst < 1e-3, sprintf("(%d LRTs, max |dLR| %.2g)", 2 * nrow(res), worst))

lead = res[res$cis_snp %in% rownames(G)[(0:(n_gene - 1)) * n_snp + 1], ]
check("G x condition detected at the causal SNPs", all(lead$p_interact < 0.05))
z = (lead$beta_geno - 0.4) / lead$se_geno                  # marginal slope at mean(cond) ~ 0
check("GLS slope consistent with the truth", all(abs(z) < 4), sprintf("(max |z| %.2f)", max(abs(z))))
check("bootstrap columns", all(c("l_boot_geno", "l_boot_interact") %in% names(res)) &&
        mean(!is.na(res$l_boot_interact)) > 0.95)

# reproducibility across core counts, and checkpoint files
run = function(cores, ...) {
  registerDoMC(cores); set.seed(11)
  r = map_interaction_qtl(Y[1:4, ], G, geneloc, snploc, anno, sample_kernel = K, normalization_approach = "none",
                          permutation_approach = "boot", cisdist = 100, ...)
  registerDoMC(1); as.data.frame(r)
}
r1 = run(1); r2 = run(2)
check("same set.seed => same bootstrap draws on 1 and 2 cores", identical(r1, r2))
cdir = file.path(tempdir(), "suez_ckpt")                   # no trailing slash
c1 = run(1, checkpoint_dir = cdir)
check("checkpoint files written inside checkpoint_dir", setequal(dir(cdir), paste0("g", 1:4, ".txt.gz")))
c2 = run(1, checkpoint_dir = cdir)                          # read back from the files
num = sapply(c1, is.numeric)
check("checkpointed results read back", identical(c1$cis_snp, c2$cis_snp) &&
        max(abs(as.matrix(c1[num]) - as.matrix(c2[num])), na.rm = TRUE) < 1e-8)
unlink(cdir, recursive = TRUE)

# maf_threshold filters on the minor allele frequency
Gm = G; s_hi = rownames(G)[2]; s_lo = rownames(G)[3]
Gm[s_hi, ] = 2; Gm[s_hi, 1:4] = 1                          # ALT frequency ~0.98: minor allele ~0.02
Gm[s_lo, ] = 0; Gm[s_lo, 1:8] = 1                          # ALT frequency ~0.03 (8 alleles: not the < 5 skip)
rm_ = map_interaction_qtl(Y[1, , drop = FALSE], Gm, geneloc, snploc, anno, sample_kernel = K,
                          normalization_approach = "none", permutation_approach = "none", cisdist = 100,
                          maf_threshold = 0.05)
check("maf_threshold drops low-MAF SNPs coded either way", !any(c(s_hi, s_lo) %in% rm_$cis_snp) &&
        all(setdiff(rownames(G)[1:n_snp], c(s_hi, s_lo)) %in% rm_$cis_snp))

# errors
check("non-PSD kernel refused", grepl("not positive semi-definite",
      errmsg(map_interaction_qtl(Y, G, geneloc, snploc, anno, sample_kernel = K - 3 * diag(N),
                                 normalization_approach = "none", permutation_approach = "none", cisdist = 100))))
one_rep = c(inds, inds[1])                                 # one repeated individual: one within contrast
a1 = data.frame(individual = one_rep, condition = rnorm(length(one_rep)))
Y1 = matrix(rnorm(2 * length(one_rep)), 2, dimnames = list(c("g1", "g2"), paste0("t", seq_along(one_rep))))
check("unidentifiable residual variance refused", grepl("not identifiable",
      errmsg(map_interaction_qtl(Y1, G, geneloc, snploc, a1, normalization_approach = "none",
                                 permutation_approach = "none", cisdist = 100))))

cat(fails, "failure(s)\n")
if (fails > 0) stop(fails, " check(s) failed")
