# Runs every usage example (doc/src/examples.R) and exports what the usage page plots: doc/src/usage_data.json.
#   Rscript doc/src/usage_data.R        (about a minute; then python3 doc/src/build.py)
suppressPackageStartupMessages(library(jsonlite))
src = if (file.exists("doc/src/examples.R")) "doc/src" else dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
set.seed(11)
source(file.path(src, "examples.R"), echo = FALSE)
registerDoMC(4)
r3 = function(x) round(x, 3)
nl = function(p) r3(-log10(pmax(p, 1e-300)))
qq_pts = function(p, n = 400) {   # expected / observed -log10 p, thinned to n quantiles plus the 30 largest
  p = sort(p[!is.na(p)]); m = length(p); e = -log10(ppoints(m)); o = -log10(p)
  i = sort(unique(c(round(seq(1, m, length.out = n)), 1:min(30, m))))
  list(e = r3(e[i]), o = r3(o[i]), lambda = r3(lambda(p)))
}
truth = read.delim(ext("truth.tsv"))
out = list()

# step 1: variance components, latent-factor variances, and a block of the learned covariance
pf = step1$optimized_fit$par
out$step1 = list(s = as.list(r3(pf$s)), factors = r3(sort(pf$s_adjustable, decreasing = TRUE)),
                 cov = r3(step1$cov_suez[1:16, 1:16]), labels = samples$sample_id[1:16])

# two conditions: top response eQTL, expression by genotype and condition, regional plot
top = regenes[1, ]
g = dosage[top$snp, samples$individual]
d = res[res$gene == top$gene, ]
out$two = list(gene = top$gene, snp = top$snp, planted = truth$type[match(top$gene, truth$gene_id)],
               expr = r3(logcpm[top$gene, ]), geno = unname(g), cond = samples$condition,
               p_geno = signif(res$p_geno[res$gene == top$gene & res$cis_snp == top$snp], 3), p_interact = signif(top$p, 3),
               pos = snploc$pos[match(d$cis_snp, snploc$snpid)], lg = nl(d$p_geno), li = nl(d$p_interact),
               body = unname(unlist(geneloc[geneloc$geneid == top$gene, c("left", "right")])), lead = top$snp, snps = d$cis_snp)

# eQTL only: top eQTL in the control samples
tc = res_ctrl[which.min(res_ctrl$p_geno), ]
out$eqtl = list(gene = tc$gene, snp = tc$cis_snp, expr = r3(logcpm[tc$gene, ctrl]), geno = unname(dosage[tc$cis_snp, anno_ctrl$individual]),
                p = signif(tc$p_geno, 3), beta = r3(tc$beta_geno), se = r3(tc$se_geno),
                planted = truth$type[match(tc$gene, truth$gene_id)])

# simulated study: genes G1..G40, causal SNP = first of each gene; effect type by gene index
type_of = function(gene) c("constant eQTL", "6h/24h response", "dose response", "mid-life peak")[as.integer(sub("G", "", gene)) %% 4 + 1]
causal = rownames(geno)[(0:(n_gene - 1)) * n_snp + 1]
lead_of = function(r, type, col = "p_interact") { r = r[r$cis_snp %in% causal & type_of(r$gene) == type, ]; r[which.min(r[[col]]), ] }

ft = lead_of(res_tp, "6h/24h response")
gt = geno[ft$cis_snp, meta$individual]
out$factor = list(gene = ft$gene, snp = ft$cis_snp, p = signif(ft$p_interact, 3), df = ft$df_interact,
                  expr = r3(expr[ft$gene, ]), geno = unname(gt), tp = as.character(meta$timepoint))

fd = lead_of(res_dose, "dose response")
out$numeric = list(gene = fd$gene, snp = fd$cis_snp, p = signif(fd$p_interact, 3), x = r3(anno_dose$condition),
                   expr = r3(expr[fd$gene, ]), geno = unname(geno[fd$cis_snp, meta$individual]))

fa = lead_of(res_age, "mid-life peak")
anno_lin = data.frame(individual = meta$individual, condition = (meta$age - 55) / 20)
rl = map_interaction_qtl(expr[fa$gene, , drop = FALSE], geno, geneloc_t, snploc_t, anno_lin,
                         curve_at = matrix((ages - 55) / 20, dimnames = list(ages, NULL)), permutation_approach = "none")
rl = rl[rl$cis_snp == fa$cis_snp, ]
out$spline = list(gene = fa$gene, snp = fa$cis_snp, ages = ages, p_spline = signif(fa$p_interact, 3), p_linear = signif(rl$p_interact, 3),
                  beta = r3(unlist(fa[paste0("beta_curve_", ages)])), se = r3(unlist(fa[paste0("se_beta_curve_", ages)])),
                  beta_lin = r3(unlist(rl[paste0("beta_curve_", ages)])), truth = r3(0.7 * exp(-((ages - 55) / 12)^2)),
                  mu = r3(unlist(fa[paste0("mu_curve_", ages)])), se_mu = r3(unlist(fa[paste0("se_mu_curve_", ages)])))

m = merge(res_tp[, c("gene", "cis_snp", "p_interact")], res_u[, c("gene", "cis_snp", "p_interact")], by = c("gene", "cis_snp"))
out$unbalanced = list(full = nl(m$p_interact.x), sub = nl(m$p_interact.y), causal = m$cis_snp %in% causal & type_of(m$gene) == "6h/24h response",
                      per_ind = as.list(table(table(meta$individual[keep]))), n_full = nrow(meta), n_sub = length(keep))

mc = merge(res_tp[, c("gene", "cis_snp", "p_geno")], res_cov[, c("gene", "cis_snp", "p_geno")], by = c("gene", "cis_snp"))
out$covariates = list(without = nl(mc$p_geno.x), with = nl(mc$p_geno.y), causal = mc$cis_snp %in% causal,
                      batch_means = as.list(r3(tapply(colMeans(expr), meta$batch, mean))))

out$related = list(cov = r3(cov_tp$cov_suez[1:18, 1:18]), grm = r3(K[1:18, 1:18]), labels = colnames(expr)[1:18])

# local ancestry: per-ancestry effects at each gene's lead SNP, and planted EUR eQTLs found pooled vs by ancestry
truth_a = read.delim(adm("truth.tsv"))
lead = do.call(rbind, lapply(split(res_a, res_a$gene), function(d) d[which.min(d$p_geno), ]))
lead = lead[lead$p_geno < 1e-3 & !is.na(lead$beta_AFR) & !is.na(lead$beta_EUR), ]
at = merge(truth_a, res_a, by.x = c("gene_id", "snp"), by.y = c("gene", "cis_snp"))
ap = merge(truth_a, res_pooled, by.x = c("gene_id", "snp"), by.y = c("gene", "cis_snp"))
found = sapply(split(at$p_geno < 1e-3, at$type), mean); found_p = sapply(split(ap$p_geno < 1e-3, ap$type), mean)
out$la = list(b_afr = r3(lead$beta_AFR), b_eur = r3(lead$beta_EUR), type = ifelse(is.na(m0 <- truth_a$type[match(lead$gene, truth_a$gene_id)]), "none planted", m0),
              found = as.list(r3(found)), found_pooled = as.list(r3(found_p)))

out$fdr = list(geno = nl(egenes$p_gene), geno_sig = egenes$fdr < 0.05, inter = nl(regenes$p_gene), inter_sig = regenes$fdr < 0.05,
               planted_inter = !is.na(truth$type[match(regenes$gene, truth$gene_id)]) & truth$type[match(regenes$gene, truth$gene_id)] == "response eQTL")
out$calib = list(geno = qq_pts(res$p_geno), inter = qq_pts(res$p_interact), boot = qq_pts(res$p_interact_boot))
out$perm = list(geno = qq_pts(res_null$p_geno), inter = qq_pts(res_null$p_interact))

writeLines(toJSON(out, auto_unbox = TRUE, digits = NA), file.path(src, "usage_data.json"))
cat("wrote", file.path(src, "usage_data.json"), "\n")
