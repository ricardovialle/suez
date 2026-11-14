# per-segment R cost of the gene loop body (base12 code, P_many, 1000 genes), 1 core
suppressPackageStartupMessages({ library(suez, lib.loc = "bench/lib/base12"); library(dplyr) })
ns = asNamespace("suez"); d = readRDS("bench/data/P_many.rds"); g = head(d$geneloc$geneid, 1000)
input = ns$quantile_normalize(d$input[g, ]); geneloc = as.data.frame(d$geneloc); snploc = d$snploc; genotype = d$genotype; anno = d$anno
eig = ns$kernel_contrast_eigen(tcrossprod(outer(anno$individual, unique(anno$individual), "==") * 1)); U = eig$vectors; lambda = eig$values
C = model.matrix(~ condition, data = anno)[, -1, drop = FALSE]; X0 = ns$standardize_cols(C); UX0 = crossprod(U, X0)
geno_row = match(as.character(snploc$snpid), rownames(genotype))
by_chr = lapply(split(seq_len(nrow(snploc)), as.character(snploc$chr)), function(i) i[order(snploc$pos[i])])
cis_rows = function(chrom, left, right) { i = by_chr[[as.character(chrom)]]; if (is.null(i)) return(integer(0))
  p = snploc$pos[i]; from = findInterval(left, p)
  r = sort(i[seq.int(from + 1L, length.out = max(0L, findInterval(right, p, left.open = TRUE) - from))]); geno_row[r[!is.na(geno_row[r])]] }
rownames(geneloc) = geneloc$geneid; cisdist = 1e5
T = numeric(0); tic = function(k) { t1 = proc.time()[[3]]; T[k] <<- (if (is.na(T[k])) 0 else T[k]) + t1 - t0; t0 <<- t1 }
res = vector("list", length(g)); set.seed(1)
for (j in seq_along(g)) { gene = g[j]; t0 = proc.time()[[3]]
  rows = cis_rows(geneloc[gene,"chr"], geneloc[gene,"left"]-cisdist, geneloc[gene,"right"]+cisdist); tic("cis_rows+geneloc[]")
  geno_here = genotype[rows,,drop=F]; tic("geno subset")
  imp_geno = if (anyNA(geno_here)) ns$easy_impute(geno_here) else { g2 = round(geno_here); storage.mode(g2) = "integer"; g2 }; tic("round/impute")
  cat(""); tic("cat")
  G = t(imp_geno[, anno$individual, drop=F]); storage.mode(G) = "double"; tic("t(G[,ind])")
  keep = rowSums(imp_geno) >= 5 & colSums(G != rep(G[1, ], each = nrow(G))) > 0; G = G[, keep, drop=F]; tic("keep")
  Uy = drop(crossprod(U, as.numeric(input[gene,]))); tic("Uy")
  fit = .Call(ns$suez_fit_gene, U, lambda, Uy, UX0, G, C, TRUE); tic(".Call fit")
  m = fit$res
  gr = data.frame(gene=gene, cis_snp=colnames(G), l0=fit$l0, l_geno=m[,1], l_interact=m[,2], stringsAsFactors=F)
  gr$p_geno = pchisq(pmax(0, 2*(gr$l_geno-gr$l0)), 1, lower.tail=F); gr$df_interact = as.integer(m[,3])
  gr$p_interact = pchisq(pmax(0, 2*(gr$l_interact-gr$l_geno)), m[,3], lower.tail=F)
  gr$beta_geno = m[,4]; gr$se_geno = m[,5]; gr$l_boot_geno = m[,6]; gr$l_boot_interact = m[,7]; tic("data.frame+pchisq")
  res[[j]] = gr }
t0 = proc.time()[[3]]; out = bind_rows(res); tic("bind_rows (once)")
print(round(cbind(total_s = T, us_per_gene = 1e6 * T / length(g)), 4))
