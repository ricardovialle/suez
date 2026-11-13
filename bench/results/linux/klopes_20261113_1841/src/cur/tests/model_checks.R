# Checks of map_interaction_qtl's model fits (run by R CMD check, or Rscript tests/model_checks.R):
#   1. every likelihood-ratio test equals the exact profile-ML LRT (grid + optimize over
#      h = s2k / (s2k + s2), b and the scale profiled out) on the intercept-free contrasts, with
#      repeated samples, a centred GRM (K 1 = 0) and a condition varying within individuals,
#      including genes whose ML sits on the s2 = 0 or s2k = 0 boundary
#   2. GLS effect sizes, bootstrap columns
#   3. a bimodal likelihood in h (two near-equal peaks)
#   4. errors: non-PSD kernel, unidentifiable residual variance
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

# missing genotypes, including a SNP constant among its calls (NaN after scaling): imputed, the gene is kept
Gi = G; Gi[cbind(sample(nrow(G), 30, TRUE), sample(nd, 30, TRUE))] = NA
s_mono = rownames(G)[2]; Gi[s_mono, ] = 0L; Gi[s_mono, 1:3] = NA
ri = map_interaction_qtl(Y[1, , drop = FALSE], Gi, geneloc, snploc, anno, sample_kernel = K, normalization_approach = "none",
                         permutation_approach = "none", cisdist = 100)
check("missing genotypes imputed, gene with a monomorphic SNP with missing calls kept",
      nrow(ri) == n_snp - 1 && !s_mono %in% ri$cis_snp && !anyNA(ri$l_geno))

# bimodal likelihood in h: two peaks of near-equal height, the higher one on either side, so a fit
# that refines only the best grid point lands on the wrong peak in one of the two genes
exact_ml = function(UX, Uy, lam) {                         # dense log-delta grid, every local max refined
  f = function(x) { h = 1 / (1 + 10^-x); d = h * lam + (1 - h); W = 1 / d
    b = solve(crossprod(UX, UX * W), crossprod(UX, Uy * W))
    -0.5 * length(Uy) - 0.5 * sum(log(sum((Uy - UX %*% b)^2 * W) / length(Uy) * d)) }
  xs = seq(-8, 8, length.out = 1601); v = sapply(xs, f); pk = which(diff(sign(diff(v))) == -2) + 1
  sapply(pk, function(i) optimize(f, xs[c(i - 1, i + 1)], maximum = TRUE, tol = 1e-10)$objective)
}
set.seed(1859)
nb = 120; Qb = qr.Q(qr(matrix(1, nb, 1)), complete = TRUE)[, -1]
lb = c(rep(1e-3, 40), rep(1, 40), rep(30, 39)) * exp(rnorm(nb - 1, 0, 0.3))
vb = c(runif(1, 0.01, 2), runif(1, 0.01, 2), runif(1, 0.01, 60))
zb = rnorm(nb - 1) * sqrt(rep(vb, c(40, 40, 39)))
Kb = Qb %*% (lb * t(Qb)) + 1 / nb                           # contrast eigenvalues lb
eb = suez:::kernel_contrast_eigen(Kb); condb = rnorm(nb); UX0b = crossprod(eb$vectors, std(condb))
y_of = function(t) drop(Qb %*% (zb * rep(c(1, 1, t), c(40, 40, 39))))
peaks = function(t) exact_ml(UX0b, drop(crossprod(eb$vectors, y_of(t))), eb$values)
sgap = function(t) { p = peaks(t); p[1] - p[length(p)] }   # left peak minus right peak
tb = sapply(c(2e-3, -2e-3), function(g) uniroot(function(t) sgap(t) - g, c(0.9, 1))$root)
indb = paste0("B", 1:nb); sb = paste0("b", 1:nb)
Yb = rbind(bL = y_of(tb[1]), bR = y_of(tb[2])); colnames(Yb) = sb; dimnames(Kb) = list(sb, sb)
Gb = matrix(rbinom(4 * nb, 2, 0.3), 4, nb, dimnames = list(paste0("sb", 1:4), indb))
rb = map_interaction_qtl(Yb, Gb, data.frame(geneid = rownames(Yb), chr = "1", left = 100, right = 110),
                         data.frame(snpid = rownames(Gb), chr = "1", pos = 101:104),
                         data.frame(individual = indb, condition = condb), sample_kernel = Kb,
                         normalization_approach = "none", permutation_approach = "none", cisdist = 100)
db = sapply(seq_len(nrow(rb)), function(i) {
  Uy = drop(crossprod(eb$vectors, Yb[rb$gene[i], ]))
  p0 = exact_ml(UX0b, Uy, eb$values)
  p1 = exact_ml(cbind(UX0b, crossprod(eb$vectors, std(Gb[rb$cis_snp[i], ]))), Uy, eb$values)
  c(length(p0), rb$l0[i] - max(p0), rb$l_geno[i] - max(p1))
})
check("bimodal likelihood: the higher of two near-equal peaks is found",
      all(db[1, ] == 2) && max(abs(db[2:3, ])) < 1e-6, sprintf("(max |dl| %.2g)", max(abs(db[2:3, ]))))

# kernel scale: the likelihood is invariant (only delta = s2k / s2 moves), including optima pushed beyond the
# grid's |log10 delta| <= 5 and zero eigenvalues with near-zero rounding noise (h -> 1)
sc = lapply(c(1, 1e6, 1e-6), function(f) map_interaction_qtl(Y, G, geneloc, snploc, anno, sample_kernel = f * K,
  normalization_approach = "none", permutation_approach = "none", cisdist = 100))
dsc = max(sapply(sc[-1], function(r) max(abs(as.matrix(r[c("l0", "l_geno", "l_interact")]) - as.matrix(sc[[1]][c("l0", "l_geno", "l_interact")])))))
check("kernel scaled by 1e6 and 1e-6: same log-likelihoods", dsc < 1e-6, sprintf("(max |dl| %.2g)", dsc))

# eigen solver (options(suez.eigen_solver)): Eigen's and LAPACK's bases differ only inside tied eigenvalues
es = lapply(c("eigen", "lapack"), function(s) { op = options(suez.eigen_solver = s); on.exit(options(op))
  map_interaction_qtl(Y, G, geneloc, snploc, anno, sample_kernel = K, normalization_approach = "none",
                      permutation_approach = "none", cisdist = 100) })
des = max(abs(as.matrix(es[[1]][c("l0", "l_geno", "l_interact")]) - as.matrix(es[[2]][c("l0", "l_geno", "l_interact")])))
check("Eigen and LAPACK eigen solvers: same log-likelihoods", des < 1e-8, sprintf("(max |dl| %.2g)", des))

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

e1 = suez:::kernel_contrast_eigen(K); G1 = matrix(rbinom(N * 3, 2, 0.3), N) * 1.0; G1[2, 2] = NaN
check("non-finite genotypes refused (no sort on NaN)", grepl("non-finite", errmsg(.Call(suez:::suez_fit_gene,
      e1$vectors, e1$values, drop(crossprod(e1$vectors, rnorm(N))), crossprod(e1$vectors, matrix(1, N, 1)),
      G1, matrix(0, N, 0), FALSE))))

cat(fails, "failure(s)\n")
if (fails > 0) stop(fails, " check(s) failed")
