# Statistical validation of a candidate suez build against the baseline.
#
# CLI (run from the repo root):
#   Rscript bench/proto/review/validate.R compare <a.rds> <b.rds> [dataset]   # a = reference (baseline)
#   Rscript bench/proto/review/validate.R truth   <res.rds> <dataset>         # power / FDR vs d$truth
#   Rscript bench/proto/review/validate.R null    <lib> <dataset> [n_genes] [out.rds] [seed]
#        runs map_interaction_qtl of bench/lib/<lib> on the dataset with individuals' genotypes permuted
#        (all associations broken, LD kept) and reports p-value calibration.
#   Rscript bench/proto/review/validate.R calib   <res.rds>                   # calibration of a null run
# From R: source("bench/proto/review/validate.R"); validate(a, b, truth = d$truth)
#   a, b: result tables (data.frame or path to RDS from bench/bench.sh / bench/run.R), or harness_table()
#   of a bench/proto/harness.R run_impl() output, so kernel-level prototypes are validated the same way.

TOL = list(abs_ll = 1e-6, rel_beta = 1e-6)   # CONTEXT.md tolerances
LL_COLS = c("l0", "l_geno", "l_interact", "l_boot_geno", "l_boot_interact", "l_boot_interact_null")

read_tab = function(x) { if (is.character(x)) x = readRDS(x); as.data.frame(x, stringsAsFactors = FALSE) }

# harness run_impl() output (+ the prepare() list) -> table with map_interaction_qtl's columns
harness_table = function(run, P) {
  do.call(rbind, lapply(seq_along(run$out), function(i) {
    f = run$out[[i]]; m = f$res; g = P$genes[[i]]
    r = data.frame(gene = g$gene, cis_snp = if (is.null(colnames(g$G))) seq_len(ncol(g$G)) else colnames(g$G),
                   l0 = f$l0, l_geno = m[, 1], l_interact = m[, 2], stringsAsFactors = FALSE)
    r$p_geno = pchisq(pmax(0, 2 * (r$l_geno - r$l0)), 1, lower.tail = FALSE)
    r$df_interact = as.integer(m[, 3])
    r$p_interact = pchisq(pmax(0, 2 * (r$l_interact - r$l_geno)), m[, 3], lower.tail = FALSE)
    r$beta_geno = m[, 4]; r$se_geno = m[, 5]; r$l_boot_geno = m[, 6]; r$l_boot_interact = m[, 7]
    r
  }))
}

# README's pooled bootstrap p-value for the interaction test
add_boot_p = function(r) {
  if (!all(c("l_boot_geno", "l_boot_interact", "l_interact") %in% names(r))) return(r)
  lr = 2 * (r$l_interact - r$l_geno); ln = sort(2 * (r$l_boot_interact - r$l_boot_geno))
  r$p_interact_boot = ifelse(is.na(lr), NA, (length(ln) - findInterval(lr, ln, left.open = TRUE) + 1) / (length(ln) + 1))
  r
}

lead_p = function(r, pc) {   # per gene: lead SNP (min p, first on ties), its p, number of SNPs tested
  r = r[!is.na(r[[pc]]), ]
  o = r[order(r$gene, r[[pc]]), ]; o = o[!duplicated(o$gene), ]
  data.frame(gene = o$gene, snp = o$cis_snp, p = o[[pc]], m = as.vector(table(r$gene)[o$gene]), stringsAsFactors = FALSE)
}
egenes = function(r, pc, fdr = 0.05) {   # Bonferroni within gene, BH across genes
  l = lead_p(r, pc); l$gene[p.adjust(pmin(1, l$p * l$m), "BH") < fdr]
}

validate = function(a, b, truth = NULL, label = "", thresholds = c(1e-3, 1e-5, 1e-8), quiet = FALSE) {
  a = add_boot_p(read_tab(a)); b = add_boot_p(read_tab(b))
  ka = paste(a$gene, a$cis_snp); kb = paste(b$gene, b$cis_snp)
  rows = c(only_a = sum(!ka %in% kb), only_b = sum(!kb %in% ka), dup_a = sum(duplicated(ka)), dup_b = sum(duplicated(kb)))
  b = b[match(ka, kb), ]; keep = !is.na(b$gene); a = a[keep, ]; b = b[keep, ]
  num = intersect(names(a)[sapply(a, is.numeric)], names(b)[sapply(b, is.numeric)])
  cols = do.call(rbind, lapply(num, function(cn) {
    x = a[[cn]]; y = b[[cn]]; ok = !is.na(x) & !is.na(y)
    d = abs(x - y)[ok]; den = pmax(abs(x), abs(y))[ok]
    if (cn == "beta_geno" && !is.null(a$se_geno)) den = pmax(den, a$se_geno[ok], na.rm = TRUE)   # beta ~ 0: relative to its se
    rel = ifelse(d == 0, 0, d / den)
    lp = if (startsWith(cn, "p_")) max(c(0, abs(log10(x[ok]) - log10(y[ok]))[x[ok] > 0 & y[ok] > 0])) else NA
    tol_fail = if (cn %in% LL_COLS) sum(d > TOL$abs_ll) else if (cn %in% c("beta_geno", "se_geno")) sum(rel > TOL$rel_beta) else
      if (cn == "df_interact") sum(d != 0) else NA
    data.frame(column = cn, n = sum(ok), na_mismatch = sum(is.na(x) != is.na(y)), identical = sum(d == 0),
               max_abs = if (any(ok)) max(d) else NA, max_rel = if (any(ok)) max(rel) else NA,
               max_dlog10p = lp, n_over_tol = tol_fail)
  }))
  pcs = intersect(c("p_geno", "p_interact", "p_interact_boot", "p_het", "p_interact_het"), num)
  genes = unique(a$gene)
  rank = do.call(rbind, lapply(pcs, function(pc) {
    la = lead_p(a, pc); lb = lead_p(b, pc); lb = lb[match(la$gene, lb$gene), ]
    sp = sapply(split(seq_len(nrow(a)), a$gene), function(i) {
      x = a[[pc]][i]; y = b[[pc]][i]; ok = !is.na(x) & !is.na(y)
      if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) NA else cor(x[ok], y[ok], method = "spearman") })
    top5 = sapply(split(seq_len(nrow(a)), a$gene), function(i) {
      setequal(head(a$cis_snp[i][order(a[[pc]][i])], 5), head(b$cis_snp[i][order(b[[pc]][i])], 5)) })
    ea = egenes(a, pc); eb = egenes(b, pc)
    sig = sapply(thresholds, function(t) { x = which(a[[pc]] < t); y = which(b[[pc]] < t)
      c(a = length(x), b = length(y), both = length(intersect(x, y))) })
    data.frame(p = pc, genes = length(la$gene), lead_same = sum(la$snp == lb$snp, na.rm = TRUE),
               top5_same = sum(top5), spearman_min = suppressWarnings(min(sp, na.rm = TRUE)),
               spearman_median = median(sp, na.rm = TRUE),
               sig = paste(sprintf("p<%g:%d/%d/%d", thresholds, sig[1, ], sig[2, ], sig[3, ]), collapse = " "),
               egenes_fdr5 = sprintf("%d/%d/%d", length(ea), length(eb), length(intersect(ea, eb))))
  }))
  equivalent = all(rows == 0) && all(cols$na_mismatch == 0) && all(cols$n_over_tol == 0, na.rm = TRUE)
  out = list(label = label, rows = rows, columns = cols, ranking = rank, equivalent = equivalent,
             power_a = if (!is.null(truth)) power_fdr(a, truth), power_b = if (!is.null(truth)) power_fdr(b, truth))
  if (!quiet) print_validation(out)
  invisible(out)
}

print_validation = function(v) {
  op = options(width = 200); on.exit(options(op))
  cat("== validate", v$label, "==\nrows:", paste(names(v$rows), v$rows, sep = "=", collapse = " "), "\n")
  print(format(v$columns, digits = 3), row.names = FALSE)
  cat("ranking (sig counts a/b/both; eGenes = Bonferroni within gene + BH 5% across genes, a/b/both):\n")
  print(v$ranking, row.names = FALSE)
  if (!is.null(v$power_a)) { cat("power/FDR a:\n"); print(v$power_a, row.names = FALSE); cat("power/FDR b:\n"); print(v$power_b, row.names = FALSE) }
  cat("VERDICT:", if (v$equivalent) "EQUIVALENT within CONTEXT.md tolerances" else "NOT equivalent (see na_mismatch / n_over_tol)", "\n")
}

# gene-level power and FDR (eGenes at BH 5%; nominal lead p < 1e-5) against sim.R's truth: every causal gene
# has a genotype effect (response effects 0.4 g trt have a main effect too); only "response" genes an interaction
power_fdr = function(r, truth) {
  r = add_boot_p(r)
  do.call(rbind, lapply(intersect(c("p_geno", "p_interact", "p_interact_boot"), names(r)), function(pc) {
    pos = if (pc == "p_geno") truth$gene else truth$gene[truth$type == "response"]
    l = lead_p(r, pc)
    rbind_rule = function(rule, hit) data.frame(test = pc, rule = rule, discoveries = length(hit),
      power = mean(pos %in% hit), fdr = if (length(hit)) mean(!hit %in% pos) else 0,
      lead_is_causal = mean(l$snp[l$gene %in% intersect(hit, pos)] %in% truth$snp))
    rbind(rbind_rule("eGene BH 5%", egenes(r, pc)), rbind_rule("lead p<1e-5", l$gene[l$p < 1e-5]))
  }))
}

# calibration of a null run: all SNPs (correlated within gene: LD) and one random SNP per gene (independent)
calib = function(r, seed = 1) {
  r = read_tab(r)
  one = { set.seed(seed); i = split(seq_len(nrow(r)), r$gene); sapply(i, function(j) j[sample.int(length(j), 1)]) }
  st = function(p, what) { p = p[!is.na(p)]; n = length(p)
    data.frame(test = what, n = n, ks_p = suppressWarnings(ks.test(p, "punif")$p.value),
               lambda_gc = median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1),
               `p<0.05` = mean(p < 0.05), `p<0.01` = mean(p < 0.01), `p<1e-3` = mean(p < 1e-3),
               min_p = min(p), check.names = FALSE) }
  out = list()
  for (pc in intersect(c("p_geno", "p_interact"), names(r))) {
    out[[length(out) + 1]] = st(r[[pc]], paste(pc, "all SNPs"))
    out[[length(out) + 1]] = st(r[[pc]][one], paste(pc, "1 SNP/gene"))
  }
  if (all(c("l_boot_geno", "l_boot_interact") %in% names(r))) {   # bootstrap LR vs chi2(df) and vs the observed null LR
    df = r$df_interact; lb = 2 * (r$l_boot_interact - r$l_boot_geno)
    out[[length(out) + 1]] = st(pchisq(pmax(0, lb), df, lower.tail = FALSE), "boot LR vs chi2(df)")
    lo = 2 * (r$l_interact - r$l_geno)
    cat(sprintf("observed vs bootstrap interaction LR (null run): two-sample KS p = %.3g, mean %.3f vs %.3f\n",
                suppressWarnings(ks.test(lo, lb)$p.value), mean(lo, na.rm = TRUE), mean(lb, na.rm = TRUE)))
  }
  out = do.call(rbind, out); print(format(out, digits = 3), row.names = FALSE)
  cat("(KS p of 'all SNPs' rows is anti-conservative: SNPs in LD are not independent; judge on '1 SNP/gene')\n")
  invisible(out)
}

null_run = function(lib, dsn, n_genes = Inf, seed = 1) {
  suppressPackageStartupMessages({ library(suez, lib.loc = lib); library(doMC) }); registerDoMC(1)
  d = readRDS(file.path("bench/data", paste0(dsn, ".rds")))
  set.seed(seed); colnames(d$genotype) = colnames(d$genotype)[sample(ncol(d$genotype))]  # breaks every association
  genes = head(rownames(d$input), n_genes)
  args = list(d$input[genes, , drop = FALSE], d$genotype, d$geneloc, d$snploc, d$anno, sample_kernel = d$sample_kernel,
              normalization_approach = "qq", permutation_approach = "boot", cisdist = 1e5)
  if (!is.null(d$covariates)) args$covariates = d$covariates
  set.seed(42); as.data.frame(do.call(map_interaction_qtl, args))
}

if (sys.nframe() == 0) {
  a = commandArgs(TRUE)
  ds = function(n) readRDS(file.path("bench/data", paste0(n, ".rds")))
  switch(a[1],
    compare = validate(a[2], a[3], truth = if (length(a) >= 4) ds(a[4])$truth, label = paste(basename(a[2]), "vs", basename(a[3]))),
    truth = print(power_fdr(read_tab(a[2]), ds(a[3])$truth), row.names = FALSE),
    calib = calib(a[2]),
    null = { lib = if (grepl("/", a[2])) a[2] else file.path("bench/lib", a[2])
             r = null_run(lib, a[3], if (length(a) >= 4) as.numeric(a[4]) else Inf, if (length(a) >= 6) as.integer(a[6]) else 1)
             if (length(a) >= 5) saveRDS(r, a[5]); calib(r) },
    stop("usage: see the header of bench/proto/review/validate.R"))
}
