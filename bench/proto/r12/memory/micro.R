d = readRDS("bench/data/P_many.rds"); anno = d$anno
set.seed(1); gh = d$genotype[sample(nrow(d$genotype), 70), ]; N = nrow(anno)
tm = function(e, r = 2000) { e = substitute(e); f = parent.frame(); 1e6 * system.time(for (i in 1:r) eval(e, f))[[3]] / r }
imp = round(gh); storage.mode(imp) = "integer"; G = t(imp[, anno$individual, drop = F]); storage.mode(G) = "double"
k1 = function() rowSums(imp) >= 5 & colSums(G != rep(G[1, ], each = nrow(G))) > 0
k2 = function() rowSums(imp) >= 5 & colSums(G != rep(G[1, ], each = nrow(G))) > 0
k3 = function() { r = range; rowSums(imp) >= 5 & vapply(seq_len(ncol(G)), function(j) any(G[, j] != G[1, j]), NA) }
k4 = function() rowSums(imp) >= 5 & colSums(G) != G[1, ] * nrow(G)   # NOT equivalent in general, timing only
cat("keep base:", tm(k1()), "us; vapply any:", tm(k3()), "us\n")
cat("rowSums(imp):", tm(rowSums(imp)), " colSums(G!=rep):", tm(colSums(G != rep(G[1, ], each = nrow(G)))), "\n")
cat("round+int:", tm({g = round(gh); storage.mode(g) = "integer"}), " t(imp[,ind]):", tm({G = t(imp[, anno$individual, drop = F]); storage.mode(G) = "double"}), "\n")
ii = match(anno$individual, colnames(gh)); cat("t(gh[,ii]) double round:", tm({G2 = round(t(gh[, ii, drop = FALSE]))}), "\n")
geneloc = as.data.frame(d$geneloc); rownames(geneloc) = geneloc$geneid; gene = geneloc$geneid[500]
cat("geneloc[gene,'chr'] x3:", tm({geneloc[gene,"chr"]; geneloc[gene,"left"]; geneloc[gene,"right"]}), "\n")
gi = match(gene, geneloc$geneid); ch = geneloc$chr; cat("vector [gi] x3:", tm({ch[gi]; geneloc$left[gi]; geneloc$right[gi]}), "\n")
m = matrix(rnorm(70 * 7), 70); colnames(G) = rownames(gh); l0 = -100
f_base = function() { gr = data.frame(gene=gene, cis_snp=colnames(G), l0=l0, l_geno=m[,1], l_interact=m[,2], stringsAsFactors=F)
  gr$p_geno = pchisq(pmax(0, 2*(gr$l_geno-gr$l0)), 1, lower.tail=F); gr$df_interact = as.integer(m[,3])
  gr$p_interact = pchisq(pmax(0, 2*(gr$l_interact-gr$l_geno)), m[,3], lower.tail=F)
  gr$beta_geno = m[,4]; gr$se_geno = m[,5]; gr$l_boot_geno = m[,6]; gr$l_boot_interact = m[,7]; gr }
f_new = function() { S = nrow(m)
  l = list(gene = rep(gene, S), cis_snp = colnames(G), l0 = rep(l0, S), l_geno = m[, 1], l_interact = m[, 2],
    p_geno = pchisq(pmax(0, 2 * (m[, 1] - l0)), 1, lower.tail = FALSE), df_interact = as.integer(m[, 3]),
    p_interact = pchisq(pmax(0, 2 * (m[, 2] - m[, 1])), m[, 3], lower.tail = FALSE),
    beta_geno = m[, 4], se_geno = m[, 5], l_boot_geno = m[, 6], l_boot_interact = m[, 7])
  structure(l, class = "data.frame", row.names = c(NA, -S)) }
cat("df base:", tm(f_base()), "us; df new:", tm(f_new()), "us; identical:", identical(f_base(), f_new()), "\n")
L1 = replicate(3000, f_base(), simplify = FALSE); L2 = replicate(3000, f_new(), simplify = FALSE)
cat("bind_rows identical:", identical(dplyr::bind_rows(L1), dplyr::bind_rows(L2)), "\n")
