## Local-ancestry (Tractor-style) models for map_interaction_qtl(local_ancestry = ...).

#' Read Tractor's per-ancestry dosage and haplotype-count files
#'
#' Reads the files written by Tractor's `extract_tracts.py` / `extract_tracts_flare.py`
#' (`<prefix>.anc0.dosage.txt`, `<prefix>.anc0.hapcount.txt`, ..., optionally gzipped) into the
#' `local_ancestry` input of [map_interaction_qtl].
#'
#' @param prefix Path prefix of the files, e.g. `"cohort/ASW.phased"`.
#' @param ancestry_names Optional names for the ancestries, in file order (anc0, anc1, ...), e.g.
#'   `c("AFR", "EUR")`. Used as column suffixes in the results.
#' @return list: `dosage` and `hapcount` (lists of SNPs x individuals matrices, one per ancestry) and
#'   `snploc` (data.frame snpid, chr, pos). SNP ids are `ID`, or `chr:pos:ref:alt` when `ID` is missing
#'   (`.`) or duplicated.
#' @details Each file holds every individual and SNP of one chromosome (or chunk) as a dense matrix:
#'   read one chromosome at a time for large cohorts.
#' @importFrom readr read_tsv
#' @export
read_tractor = function(prefix, ancestry_names = NULL) {
  find = function(type) {
    f = Sys.glob(paste0(prefix, ".anc*.", type, ".txt*"))
    k = as.integer(sub(".*\\.anc([0-9]+)\\..*", "\\1", f))
    f[order(k)]
  }
  fd = find("dosage"); fh = find("hapcount")
  if (length(fd) < 2 || length(fd) != length(fh))
    stop("expected one dosage and one hapcount file per ancestry (at least 2) for prefix ", prefix, call. = FALSE)
  rd = function(f) {   # CHROM POS ID REF ALT, then one column per individual
    n = length(strsplit(readLines(f, n = 1), "\t")[[1]])
    read_tsv(f, col_types = paste0("cdccc", strrep("d", n - 5)), progress = FALSE)
  }
  first = rd(fd[1])
  meta = first[, 1:5]; names(meta) = c("chr", "pos", "id", "ref", "alt")
  id = meta$id
  bad = is.na(id) | id == "." | duplicated(id) | duplicated(id, fromLast = TRUE)
  id[bad] = paste(meta$chr, meta$pos, meta$ref, meta$alt, sep = ":")[bad]
  as_mat = function(d) { m = as.matrix(d[, -(1:5)]); rownames(m) = id; m }
  read_all = function(files) lapply(files, function(f) {
    d = if (f == fd[1]) first else rd(f)
    if (nrow(d) != nrow(meta) || any(d[[2]] != meta$pos)) stop(f, " does not match the SNPs of ", fd[1], call. = FALSE)
    as_mat(d)
  })
  dosage = read_all(fd); hapcount = read_all(fh)
  nm = if (is.null(ancestry_names)) paste0("anc", seq_along(fd) - 1) else ancestry_names
  names(dosage) = names(hapcount) = nm
  list(dosage = dosage, hapcount = hapcount,
       snploc = data.frame(snpid = id, chr = meta$chr, pos = meta$pos, stringsAsFactors = FALSE))
}

# Checks the local_ancestry input; returns it with the ancestry names set.
check_local_ancestry = function(la, individuals) {
  if (!is.list(la) || is.null(la$dosage) || is.null(la$hapcount))
    stop("local_ancestry must be a list with `dosage` and `hapcount` (see read_tractor()).", call. = FALSE)
  k = length(la$dosage)
  if (k < 2 || length(la$hapcount) != k) stop("local_ancestry needs >= 2 ancestries, with one dosage and one hapcount matrix each.", call. = FALSE)
  if (is.null(names(la$dosage))) names(la$dosage) = paste0("anc", seq_len(k) - 1)
  names(la$hapcount) = names(la$dosage)
  ref = la$dosage[[1]]
  for (m in c(la$dosage, la$hapcount)) {
    if (!is.matrix(m) || !identical(dim(m), dim(ref)) || !identical(rownames(m), rownames(ref)) || !identical(colnames(m), colnames(ref)))
      stop("all local_ancestry matrices must have the same SNPs (rows) and individuals (columns), in the same order.", call. = FALSE)
    if (anyNA(m)) stop("local_ancestry matrices must not contain NA.", call. = FALSE)
  }
  miss = setdiff(unique(individuals), colnames(ref))
  if (length(miss)) stop(length(miss), " individual(s) of anno are not columns of local_ancestry, e.g. ", miss[1], call. = FALSE)
  la
}

std_cols = function(M) {   # centre and scale; a constant column becomes 0 (= term absent for that SNP)
  M = sweep(M, 2, colMeans(M)); s = sqrt(colSums(M^2) / (nrow(M) - 1)); ok = s > 1e-10
  M[, ok] = sweep(M[, ok, drop = FALSE], 2, s[ok], "/"); M[, !ok] = 0; M
}

# One gene in local-ancestry mode. D, H: lists (one per ancestry) of cis SNPs x individuals matrices.
# Models (terms added to the condition/covariate design X0, per SNP):
#   M0  LA                        M1c LA + G + LA x C
#   M1  LA + G_1..G_k             M2  M1c + G_1..G_k x C
#   M1s LA + (G_1 + ... + G_k)    M2s M1c + (G_1 + ... + G_k) x C
# LA = haplotype counts of the first k-1 ancestries, G_a = ALT dosage on ancestry-a haplotypes.
la_gene = function(gene, D, H, anno, C, U, lambda, Uy, UX0, boot, maf_threshold) {
  k = length(D); anc = names(D); ind = anno$individual; q = ncol(C)
  Gt = Reduce(`+`, D)[, unique(ind), drop = FALSE]             # SNPs x individuals, total ALT dosage
  keep = rowSums(Gt) >= 5 & apply(Gt, 1, sd) > 0
  if (maf_threshold > 0) { af = rowMeans(Gt) / 2; keep = keep & pmin(af, 1 - af) > maf_threshold }
  if (!any(keep)) return(NULL)
  D = lapply(D, function(m) m[keep, , drop = FALSE]); H = lapply(H, function(m) m[keep, , drop = FALSE])
  ui = unique(ind)
  # SNPs x k: >= 5 copies of each allele on that ancestry's haplotypes (with all of them ALT, G_a = LA_a: collinear)
  present = sapply(seq_len(k), function(a) { alt = rowSums(D[[a]][, ui, drop = FALSE])
    alt >= 5 & rowSums(H[[a]][, ui, drop = FALSE]) - alt >= 5 })
  if (!is.matrix(present)) present = matrix(present, nrow = 1)
  Gs = lapply(seq_len(k), function(a) { g = t(D[[a]][, ind, drop = FALSE]); g[, !present[, a]] = 0; g })   # samples x S
  Ls = lapply(H[-k], function(m) t(m[, ind, drop = FALSE]))
  Gsum = t(Reduce(`+`, D)[, ind, drop = FALSE])

  cols = list(); idx = function(x) { n0 = length(cols); cols <<- c(cols, x); n0 + seq_along(x) - 1L }
  iLA = idx(lapply(Ls, std_cols)); iG = idx(lapply(Gs, std_cols)); iGs = idx(list(std_cols(Gsum)))
  iLAc = iGc = iGsc = integer(0)
  for (j in seq_len(q)) {
    iLAc = c(iLAc, idx(lapply(Ls, function(x) std_cols(x * C[, j]))))
    iGc = c(iGc, idx(lapply(Gs, function(x) std_cols(x * C[, j]))))
    iGsc = c(iGsc, idx(list(std_cols(Gsum * C[, j]))))
  }
  models = list(M0 = iLA, M1 = c(iLA, iG), M1s = c(iLA, iGs), M1c = c(iLA, iG, iLAc),
                M2 = c(iLA, iG, iLAc, iGc), M2s = c(iLA, iG, iLAc, iGsc))
  do_boot = boot && q > 0
  fit = .Call(suez_fit_models, U, lambda, Uy, UX0, cols, models, if (do_boot) 3L else -1L, if (do_boot) 4L else -1L)
  ll = fit$ll; np = fit$npar; colnames(ll) = colnames(np) = names(models)
  test = function(a, b) {
    df = np[, a] - np[, b]
    list(df = as.integer(df), p = ifelse(df > 0, pchisq(pmax(0, 2 * (ll[, a] - ll[, b])), pmax(df, 1), lower.tail = FALSE), NA_real_))
  }
  tg = test("M1", "M0"); th = test("M1", "M1s"); ti = test("M2", "M1c"); tih = test("M2", "M2s")
  out = data.frame(gene = gene, cis_snp = rownames(D[[1]]), l0 = ll[, "M0"], l_geno = ll[, "M1"],
                   l_geno_shared = ll[, "M1s"], l_interact_null = ll[, "M1c"], l_interact = ll[, "M2"],
                   l_interact_shared = ll[, "M2s"], df_geno = tg$df, p_geno = tg$p, df_het = th$df, p_het = th$p,
                   df_interact = ti$df, p_interact = ti$p, p_interact_het = tih$p, stringsAsFactors = FALSE)
  if (q == 0) out[c("l_interact_null", "l_interact", "l_interact_shared", "p_interact", "p_interact_het")] = NA_real_
  # Wald SEs use the residual variance on nc - p degrees of freedom and a t reference, as lm() (and
  # Tractor) do: with the default kernel and one sample per individual the results equal Tractor's.
  nc = length(Uy)
  rdf = function(m) pmax(nc - ncol(UX0) - np[, m], 1)
  wald = function(b, s, df) 2 * pt(-abs(b / s), df)
  df1 = rdf("M1"); df2 = rdf("M2")
  for (a in seq_len(k)) {
    sdr = apply(Gs[[a]], 2, sd)
    b = fit$beta[[2]][, k - 1 + a] / sdr; s = fit$se[[2]][, k - 1 + a] / sdr * sqrt(nc / df1)
    hc = rowSums(H[[a]][, ui, drop = FALSE])
    out[[paste0("af_", anc[a])]] = ifelse(hc > 0, rowSums(D[[a]][, ui, drop = FALSE]) / hc, NA_real_)
    out[[paste0("laprop_", anc[a])]] = hc / (2 * length(ui))
    out[[paste0("beta_", anc[a])]] = b; out[[paste0("se_", anc[a])]] = s; out[[paste0("p_", anc[a])]] = wald(b, s, df1)
    if (q == 1) {
      pos = (k - 1) + k + (k - 1) + a                           # G_a x C within M2
      sdc = apply(Gs[[a]] * C[, 1], 2, sd)
      bi = fit$beta[[5]][, pos] / sdc; si = fit$se[[5]][, pos] / sdc * sqrt(nc / df2)
      out[[paste0("beta_interact_", anc[a])]] = bi; out[[paste0("se_interact_", anc[a])]] = si
      out[[paste0("p_interact_", anc[a])]] = wald(bi, si, df2)
    }
  }
  if (do_boot) { out$l_boot_interact_null = fit$ll_boot[, 1]; out$l_boot_interact = fit$ll_boot[, 2] }
  out
}
