# Genotype calls (SNPs x individuals) rounded to allele counts; a missing call gets its SNP's mean count over the
# calls (mean dosage, as tensorQTL, FastQTL / QTLtools and GEMMA), the same in every gene's window. A SNP with no
# calls is set constant (0), so the variance filter drops it.
impute_mean = function(geno) {
  g = round(geno)
  if (anyNA(g)) {
    ij = which(is.na(g), arr.ind = TRUE); mu = rowMeans(g, na.rm = TRUE); mu[is.nan(mu)] = 0
    g[ij] = mu[ij[, 1]]
  }
  g
}

#' qqnorm without plotting
#'
#' @param x Numeric vector input
#' @return Quantile normalized version
qqnorm_no_plot=function(x) {
  n=length(x)
  a=if (n<=10) (3.0/8.0) else 0.5
  qnorm( (rank(x)-a)/(n+1.0-2.0*a) )
}

#' Quantile normalize columns
#' @param input Numeric matrix
quantile_normalize_cols=function(input) {
  apply(input, 2, qqnorm_no_plot)
}

#' Quantile normalize rows (to normal)
#' @param input Numeric matrix
#' @importFrom magrittr %>%
quantile_normalize=function(input) {
  # row by row into one result matrix: apply() + t() hold ~6 copies of the input at once
  out = matrix(0, nrow(input), ncol(input), dimnames = dimnames(input))
  for (i in seq_len(nrow(input))) out[i, ] = qqnorm_no_plot(input[i, ])
  out
}

get_cis_geno <- function(genotype, chrom, left, right, snploc, genome_build, ...){ UseMethod('get_cis_geno') }

#' @export
get_cis_geno.TabixFile = function(genotype, chrom, left, right, snploc, genome_build, ...) {
  for (p in c("VariantAnnotation", "GenomicRanges", "IRanges"))
    if (!requireNamespace(p, quietly = TRUE)) stop("Bioconductor package '", p, "' is needed for a TabixFile genotype.")
  gr = GenomicRanges::GRanges(chrom, IRanges::IRanges(left, right))
  sp = VariantAnnotation::ScanVcfParam(which = gr)
  vcf = VariantAnnotation::readVcf(genotype, genome_build, param = sp)
  gt = VariantAnnotation::geno(vcf)$GT
  if (nrow(gt) == 0) return(NULL)
  allele1 = substr(gt, 1, 1)
  class(allele1) = "numeric"
  allele2 = substr(gt, 3, 3)
  class(allele2) = "numeric"
  allele1 + allele2
}
                
#' @export
get_cis_geno.matrix = function(genotype, chrom, left, right, snploc, genome_build, ...) {
  cis_snps=as.character(snploc$snpid[snploc$chr == chrom & snploc$pos > left & snploc$pos < right])
  if (length(cis_snps)==0) NULL else genotype[cis_snps,,drop=F]
}

#' Response (e)QTL mapping
#'
#' @param input [genes x samples] matrix of log_2 expression values, e.g. cpm or fpkm
#' @param genotype [SNPs x individuals] matrix of genotypes (can have some NAs) OR a Bioconductor `TabixFile`` object. Ignored (can be NULL) when `local_ancestry` is given.
#' @param geneloc [genes x 4] data.frame with columns (geneid, chr, left, right). The intersect of geneloc$geneid and rownames(input) will be tested.
#' @param snploc [SNPs x 3] data.frame with columns (snpid, chr, pos) where snpid's map to rows of `genotype`
#' @param anno [samples x 2] data.frame with columns (individual, condition) where individual must correspond to column names of `genotype`. A `condition` with a single value (e.g. all 1) maps eQTLs without interaction tests.
#'   `condition` may also be a numeric matrix column, e.g. `anno$condition = splines::ns(age, df = 3)`: both the
#'   condition trajectory and the genotype effect are then smooth functions of age (varying-coefficient model),
#'   and `p_interact` tests all genotype x basis terms jointly (`df_interact` = number of columns).
#' @param sample_kernel [samples x samples] covariance matrix learned in suez step 1. If you don't want to control for latent confounders leave this as NULL to use the matrix representing which samples are from the same individual.
#' @param normalization_approach One of (qq,log,linear) determining how to normalize within each gene. We recommend qq (quantile normalization to normal) to ensure the assumptions of the test hold.
#' @param permutation_approach One of ("none","permute","boot"). We recommend boot (parametric bootstrap) since the permutation test is not really valid for interaction effect testing.
#' @param cisdist How far from the gene boundaries to look for SNPs to test
#' @param checkpoint_dir An optional directory to store per gene results. If mapping crashes these results will be reused to save time.
#' @param debug Controls whether you see progress messages and whether errors are suppressed.
#' @param genome_build Only required if `genotype` is `TabixFile`. Default is "GRCh38"
#' @param maf_threshold Variants with minor allele frequency at or below this value are not tested.
#' @param local_ancestry Optional, for admixed samples: list with `dosage` and `hapcount`, each a list
#'   with one [SNPs x individuals] matrix per ancestry (ALT dosage on that ancestry's haplotypes, and
#'   number of haplotypes of that ancestry, 0/1/2), as written by Tractor's `extract_tracts.py` and read
#'   by [read_tractor]. SNP ids (row names) map to `snploc`. Switches to the local-ancestry models below.
#' @param covariates Optional [samples x q] numeric matrix or data.frame of fixed covariates added to
#'   every model (e.g. global ancestry proportions or genotype PCs, age, sex).
#' @param curve_at Optional: condition values at which to return the fitted curves, one row per point and one
#'   column per condition term (a vector for a numeric condition; for a basis, the basis at those values, e.g.
#'   `predict(splines::ns(age, df = 3), c(30, 50, 70))`). Row names (or the values of a vector) label the
#'   columns. Adds, per point `x`: `beta_curve_x` and `se_beta_curve_x` (per-allele genotype effect at that
#'   condition, from the interaction model; in local-ancestry mode `beta_curve_<ancestry>_x`, per ancestry, SE on
#'   the residual df as the other ancestry SEs) and `mu_curve_x`, `se_mu_curve_x` (expression trajectory of the
#'   no-genotype model, the same on every row of a gene, relative to the condition whose row of `curve_at` is
#'   0, e.g. the lower boundary knot of `ns()`, since the intercept is not estimated). SEs are GLS SEs at the
#'   fitted variances.
#'
#' @details Model: y = X b + e, e ~ N(0, sigma2_k * K + sigma2 * I), K = `sample_kernel`, fitted by
#' maximum likelihood on the N-1 contrasts orthogonal to the intercept (REML for the intercept only),
#' so a kernel with K 1 = 0 (a GRM from centred genotypes) is handled. Each fit maximises the
#' likelihood exactly: b and the total variance have closed forms given h = sigma2_k / (sigma2_k +
#' sigma2), and h is found by a grid in log10(sigma2_k / sigma2) followed by Newton steps around
#' every local maximum of the grid (the likelihood can have two peaks). Failed fits give NA. Columns are standardized
#' before fitting (the likelihood is invariant). The function stops if the kernel is not positive
#' semi-definite, or if the condition design leaves no kernel null direction for sigma2 (e.g. the
#' default kernel with condition fitting every within-individual contrast).
#'
#' Parallelism: genes run on the registered doMC workers (forked processes). `options(suez.threads = n)`
#' (or the environment variable SUEZ_THREADS) also fits each gene's SNP blocks on n threads (default 1;
#' local-ancestry mode is not threaded). Results are identical at any thread count. Threads help when there are
#' few genes or N is large (e.g. a GRM kernel with N >= 1000), and use much less memory than the same number of
#' forked workers; total CPUs used = workers x threads.
#'
#' @section Local-ancestry models:
#' With `local_ancestry` (k ancestries), each SNP is fitted with the models of Tractor (Atkinson et al.
#' 2021) inside the suez mixed model: LA = haplotype counts of the first k - 1 ancestries, G_a = ALT
#' dosage on ancestry-a haplotypes, C = condition multipliers.
#' M0: LA; M1: LA + G_1..G_k; M1s: LA + (G_1 + ... + G_k); M1c: M1 + LA x C; M2: M1c + G_1..G_k x C;
#' M2s: M1c + (G_1 + ... + G_k) x C. Tests: `p_geno` (M1 vs M0, any ancestry-specific effect),
#' `p_het` (M1 vs M1s, effects differ between ancestries), `p_interact` (M2 vs M1c, any
#' ancestry-specific response; LA x C is in the null so a local-ancestry effect on the response is not
#' attributed to the SNP), `p_interact_het` (M2 vs M2s). Per ancestry: `af_*`, `laprop_*`, and GLS
#' `beta_*`, `se_*`, `p_*` (Wald, M1, per ALT allele; SE on the residual df and a t reference, as Tractor) and, for one interaction column,
#' `beta_interact_*`, `se_interact_*`, `p_interact_*` (Wald, M2); for several interaction columns (a factor
#' with more than two levels, or a basis) `p_interact_*` is the joint Wald F test of that ancestry's q
#' interaction terms (df q and the residual df). An ancestry's terms are dropped for a
#' SNP with fewer than 5 copies of either allele on that ancestry's haplotypes (its columns are NA,
#' and the df shrink).
#' The bootstrap simulates under M1c: `l_boot_interact_null`, `l_boot_interact`.
#' Differences between ancestries at a SNP that only tags the causal variant can reflect LD that
#' differs between ancestries rather than different effects.
#'
#' @return data.frame, one row per gene x SNP: `l0`, `l_geno`, `l_interact` (log-likelihoods of the
#' contrasts under the no-genotype, genotype and genotype x condition models), `p_geno` (LRT, df 1),
#' `df_interact` and `p_interact` (LRT of `l_interact` vs `l_geno`), `beta_geno` and `se_geno` (GLS
#' genotype effect per ALT allele and its SE at the fitted variances), and for
#' `permutation_approach = "boot"` `l_boot_geno` and `l_boot_interact` (one parametric bootstrap draw
#' under the genotype model: 2 * (l_boot_interact - l_boot_geno) is a null draw of the interaction LR).
#' @import doMC
#' @useDynLib suez, .registration = TRUE
#' @export
map_interaction_qtl = function(input, genotype, geneloc, snploc, anno, sample_kernel=NULL, 
                               normalization_approach="qq", permutation_approach="boot", cisdist=1e5, 
                               checkpoint_dir=NULL, debug=F, genome_build= "GRCh38", maf_threshold = 0,
                               local_ancestry=NULL, covariates=NULL, curve_at=NULL) {

  if (normalization_approach=="qq") {
    input=quantile_normalize(input)
  } else if (normalization_approach=="log") {
    input=input %>% t %>% scale %>% t
  } else if (normalization_approach=="linear") {
    input=(2^input) %>% t %>% scale %>% t
  } else if (normalization_approach !="none" ) {
    stop(paste0("Invalid normalization_approach:",normalization_approach))
  }
  
  geneloc = as.data.frame(geneloc) # code doesn't work with tibble

  cat("Checkpoint dir:", checkpoint_dir, "\n")
  if (!is.null(checkpoint_dir)) dir.create(checkpoint_dir, recursive = T, showWarnings = F)

  genes=intersect(rownames(input),geneloc$geneid)

  la_mode = !is.null(local_ancestry)
  if (la_mode) local_ancestry = check_local_ancestry(local_ancestry, anno$individual)
  ref_geno = if (la_mode) local_ancestry$dosage[[1]] else genotype
  if (is.matrix(ref_geno)) {
    # cis lookup: snploc rows per chromosome sorted by position, and snploc -> genotype rows, built once
    # (scanning all of snploc for every gene costs ~20 ms per gene per million SNPs)
    geno_row = match(as.character(snploc$snpid), rownames(ref_geno))
    by_chr = lapply(split(seq_len(nrow(snploc)), as.character(snploc$chr)), function(i) i[order(snploc$pos[i])])
    cis_rows = function(chrom, left, right) {   # genotype rows of SNPs with left < pos < right, in snploc order
      i = by_chr[[as.character(chrom)]]
      if (is.null(i)) return(integer(0))
      p = snploc$pos[i]; from = findInterval(left, p)
      r = sort(i[seq.int(from + 1L, length.out = max(0L, findInterval(right, p, left.open = TRUE) - from))])
      geno_row[r[!is.na(geno_row[r])]]
    }
  } else if (la_mode) stop("local_ancestry matrices must be matrices (see read_tractor()).", call. = FALSE)
  rownames(geneloc)=geneloc$geneid
  # per-gene inputs by position (geneloc[gene, "chr"] costs ~45 us per lookup, input[gene, ] a name match)
  gi = match(genes, geneloc$geneid)
  g_chr = geneloc$chr[gi]; g_left = geneloc$left[gi] - cisdist; g_right = geneloc$right[gi] + cisdist
  g_row = match(genes, rownames(input))



  cat("Performing initial eigendecomposition...\n")
  # N-1 contrasts orthogonal to the intercept; default kernel Z Z' (same individual): analytic sparse basis
  eig = if (is.null(sample_kernel)) zz_contrast_basis(anno$individual) else kernel_contrast_eigen(sample_kernel)
  if (is.null(sample_kernel)) eig$vectors = eig$U
  U=eig$vectors
  lambda=eig$values
  sample_kernel=NULL   # not used below: frees the default kernel (N x N) during the gene loop

  # GxE multipliers: the numeric condition (one column per column of a basis matrix, e.g. splines::ns(age)),
  # or indicators of the non-reference levels; genotype x column j is the interaction term j
  # (as in model.matrix(~ geno:condition)). A single condition (e.g. a constant) gives a plain eQTL map.
  cond_key = if (is.matrix(anno$condition)) do.call(paste, as.data.frame(anno$condition)) else as.character(anno$condition)
  one_condition = length(unique(cond_key)) == 1
  C = if (one_condition) matrix(0, nrow(anno), 0) else
    if (is.numeric(anno$condition)) cbind(anno$condition) else model.matrix(~ condition, data=anno)[, -1, drop=F]
  df_interact = ncol(C)
  sdC = apply(C, 2, sd)
  if (!is.null(curve_at)) {
    curve_at = as.matrix(curve_at)
    if (one_condition || any(sdC == 0) || !is.numeric(curve_at) || anyNA(curve_at) || ncol(curve_at) != df_interact)
      stop("curve_at must be a numeric matrix (or vector) without NA, with one column per condition term (",
           df_interact, ") of a condition that varies.", call. = FALSE)
    if (is.null(rownames(curve_at)))
      rownames(curve_at) = if (df_interact == 1) format(curve_at[, 1], trim = TRUE) else seq_len(nrow(curve_at))
  }

  no_geno =if (one_condition) matrix(0, nrow(anno), 0) else model.matrix( ~ condition, data=anno)
  if (!is.null(covariates)) {
    if (is.data.frame(covariates)) covariates = model.matrix(~ ., data = covariates)
    covariates = as.matrix(covariates)
    if (nrow(covariates) != nrow(anno) || !is.numeric(covariates) || anyNA(covariates))
      stop("covariates must be numeric, without NA, with one row per sample (as anno).", call. = FALSE)
    no_geno = cbind(no_geno, covariates)
  }
  no_geno = no_geno[, apply(no_geno, 2, sd) > 0, drop=F] # intercept: removed by the contrasts
  X0 = standardize_cols(no_geno)
  UX0 = crossprod_U(U, X0)

  # condition varying within individual => GxE columns can fit the kernel's null directions
  within = any(tapply(cond_key, anno$individual, function(x) length(unique(x)) > 1))
  k_anc = if (la_mode) length(local_ancestry$dosage) else 1
  check_identifiable(eig, X0, need = 1 + if (within) (2 * k_anc - 1) * df_interact else 0)

  threads = suez_threads()
  errorhandling=if (debug) 'stop' else 'remove'

  if (!debug) {
    zz <- file( "/dev/null", open = "wt")
    sink(zz)
    sink(zz, type = "message")
    on.exit({
      sink(type="message")
      sink()
      close(zz)
    })
  }

  # one seed per gene, drawn from the caller's RNG: results depend on set.seed() only, not on the cores
  seeds = sample.int(.Machine$integer.max, length(genes))

  # results are combined once at the end (combining pairwise copies the growing table at every gene)
  res = foreach(gene=genes, seed=seeds, j=seq_along(genes), .errorhandling=errorhandling) %dopar% {

    set.seed(seed)
    if (!is.null(checkpoint_dir)){
      check_fn = file.path(checkpoint_dir, paste0(gene, ".txt.gz"))
      if (file.exists(check_fn)) {
        return(readr::read_tsv(check_fn))
      }
    }
    
    if (la_mode) {
      rows = cis_rows(g_chr[j], g_left[j], g_right[j])
      if (!length(rows)) return(NULL)
      D = lapply(local_ancestry$dosage, function(m) m[rows,,drop=F])
      H = lapply(local_ancestry$hapcount, function(m) m[rows,,drop=F])
      if (permutation_approach=="permute") {   # same shuffle of individuals for every ancestry
        perm = sample(ncol(D[[1]]))
        D = lapply(D, function(m) { colnames(m) = colnames(m)[perm]; m }); H = lapply(H, function(m) { colnames(m) = colnames(m)[perm]; m })
      }
      Uy = drop(crossprod_U(U, as.numeric(input[g_row[j],])))
      gene_results = la_gene(gene, D, H, anno, C, U, lambda, Uy, UX0, permutation_approach=="boot", maf_threshold,
                             curve_at, sdC)
      if (!is.null(checkpoint_dir) && !is.null(gene_results)) readr::write_tsv(gene_results, check_fn)
      return(gene_results)
    }

    geno_here = if (is.matrix(genotype)) {
      rows = cis_rows(g_chr[j], g_left[j], g_right[j])
      if (length(rows)) genotype[rows,,drop=F] else NULL
    } else get_cis_geno(genotype, 
                              chrom = geneloc[gene,"chr"], 
                              left = geneloc[gene,"left"]-cisdist, 
                              right = geneloc[gene,"right"]+cisdist, 
                              snploc = snploc,
                              genome_build = genome_build)
    if (is.null(geno_here)) return(NULL)
    
    imp_geno = impute_mean(geno_here)
    if (maf_threshold > 0) {
      af = rowMeans(geno_here, na.rm=T)/2 # ALT frequency, on unimputed genotypes
      mafs = pmin(af, 1 - af)
      imp_geno = imp_geno[mafs > maf_threshold,,drop=F]
      if (nrow(imp_geno)==0) return(NULL)
    }
    
    cat(gene, nrow(imp_geno), " cis snps\n")

    if (permutation_approach=="permute") colnames(imp_geno)=colnames(imp_geno)[ sample(ncol(imp_geno),ncol(imp_geno)) ]

    G = t(imp_geno[, anno$individual, drop=F])                     # samples x SNPs
    storage.mode(G) = "double"
    # >= 5 ALT alleles; varies across samples: N sum(x^2) == (sum x)^2 iff x is constant (Cauchy-Schwarz), exact
    # here: a SNP constant among its calls is filled with that same count (integer sums divide exactly), and any
    # other column's variance is far above rounding (counts, or means of counts)
    cs = colSums(G)
    keep = rowSums(imp_geno) >= 5 & nrow(G) * colSums(G * G) != cs * cs
    if (!any(keep)) return(NULL)
    G = G[, keep, drop=F]
    Uy = drop(crossprod_U(U, as.numeric(input[g_row[j],])))

    if (debug) cat("Fitting", ncol(G), "SNPs...\n")
    curves = !is.null(curve_at)
    fit = if (threads > 1) .Call(suez_fit_gene_threads, U, lambda, Uy, UX0, G, C, permutation_approach=="boot", threads, curves)
      else .Call(suez_fit_gene, U, lambda, Uy, UX0, G, C, permutation_approach=="boot", curves)   # src/lmm.cpp
    m = fit$res
    # columns as a list (data.frame() + 8 `$<-` cost ~0.5 ms per gene)
    S = nrow(m)
    gene_results = list(gene=rep(gene, S), cis_snp=colnames(G), l0=rep(fit$l0, S), l_geno=m[,1], l_interact=m[,2],
                        p_geno=pchisq(pmax(0, 2*(m[,1]-fit$l0)), 1, lower.tail=F), df_interact=as.integer(m[,3]),
                        p_interact=pchisq(pmax(0, 2*(m[,2]-m[,1])), m[,3], lower.tail=F), beta_geno=m[,4], se_geno=m[,5])
    if (permutation_approach=="boot") {
      # parametric bootstrap under the genotype (no-interaction) model, one draw per SNP
      gene_results$l_boot_geno = m[,6]
      gene_results$l_boot_interact = m[,7]
    }
    if (curves) {
      q = ncol(C); k1 = 1 + q
      gene_results = mu_curve(gene_results, fit$b0[1:q], fit$V0[1:q, 1:q, drop = FALSE], sdC, curve_at, S)
      gene_results = add_curve(gene_results, "beta_curve_", fit$gx[, 1:k1, drop = FALSE], fit$gx[, -(1:k1), drop = FALSE],
                               cbind(1, curve_at))
    }
    gene_results = structure(gene_results, class="data.frame", row.names=c(NA, -S))

    if (!is.null(checkpoint_dir)){
      cat("Saving results to ", checkpoint_dir, "\n")
      readr::write_tsv(gene_results, check_fn)
    }

    gene_results
  }
  bind_tables(res)
}

# Curve columns <prefix><label> and se_<prefix><label>: h'b and sqrt(h'Vh) per SNP at each row h of H,
# b [S x d], V [S x d^2] column-major
add_curve = function(out, prefix, b, V, H) {
  for (p in seq_len(nrow(H))) {
    h = H[p, ]; lab = rownames(H)[p]
    out[[paste0(prefix, lab)]] = drop(b %*% h)
    out[[paste0("se_", prefix, lab)]] = sqrt(drop(V %*% as.vector(h %o% h)))
  }
  out
}

# mu_curve_<label>: condition trajectory of the no-genotype model, relative to the condition whose row of
# curve_at is 0 (the intercept is not estimated); b, V on the standardized X0 scale, repeated on S rows
mu_curve = function(out, b, V, sdC, curve_at, S) {
  add_curve(out, "mu_curve_", matrix(b / sdC, S, length(b), byrow = TRUE),
            matrix(as.vector(V / outer(sdC, sdC)), S, length(V), byrow = TRUE), curve_at)
}

# dplyr::bind_rows() for per-gene data.frames with the same columns (identical() result), without loading
# dplyr (~33 MB, ~0.15 s); anything else (checkpoint tibbles from read_tsv) goes to bind_rows
bind_tables = function(res) {
  res = res[!vapply(res, is.null, NA)]
  nm = if (length(res)) names(res[[1]])
  same = length(res) && all(vapply(res, function(x) identical(class(x), "data.frame") && identical(names(x), nm) &&
    identical(vapply(x, typeof, ""), vapply(res[[1]], typeof, "")), NA))
  if (!same) return(dplyr::bind_rows(res))
  out = lapply(seq_along(nm), function(k) unlist(lapply(res, `[[`, k), use.names = FALSE)); names(out) = nm
  structure(out, class = "data.frame", row.names = c(NA, -length(out[[1]])))
}
