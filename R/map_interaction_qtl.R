

#' Consistent implementation of diag
fix_diag=function(x) {
  if(length(x)==1) matrix(x) else diag(x)
}

unscale = function(x) {
  x = sweep(x, 1, attr(x, "scaled:scale"), "*")
  sweep(x, 1, attr(x, "scaled:center"), "+")
}

#' Simple SVD based imputation of missing genotypes
#'
#' @param geno [samples] x [SNPs] genotype matrix (0/1/2)
#' @param prop_var Proportion of variance that the PCs should explain
#'
#' @return Complete genotype matrix.
easy_impute=function(geno, prop_var=0.95) {
  temp=geno
  temp=t(scale(t(geno)))
  temp[is.na(temp)]=0
  s=svd(temp)
  v=s$d^2/sum(s$d^2)
  to_use=cumsum(v)<prop_var
  s$d[!to_use]=0.0
  recon=s$u %*% fix_diag(s$d) %*% t(s$v)
  temp[is.na(geno)]=recon[is.na(geno)]
  temp=unscale(temp)
  stopifnot(max(abs(temp[!is.na(geno)]-geno[!is.na(geno)]))<1e-10)
  temp=round(temp)
  class(temp)="integer"
  temp
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
#' `beta_interact_*`, `se_interact_*`, `p_interact_*` (Wald, M2). An ancestry's terms are dropped for a
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
                               local_ancestry=NULL, covariates=NULL) {

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

  if (is.null(sample_kernel)) {
    # same individual: Z Z' for the samples x individuals indicator Z (outer() builds two N^2 vectors of ids)
    Z = matrix(0, nrow(anno), length(unique(anno$individual)))
    Z[cbind(seq_len(nrow(anno)), match(anno$individual, unique(anno$individual)))] = 1
    sample_kernel = tcrossprod(Z)
    rm(Z)
  }

  cat("Performing initial eigendecomposition...\n")
  eig=kernel_contrast_eigen(sample_kernel) # N-1 contrasts orthogonal to the intercept
  U=eig$vectors
  lambda=eig$values
  sample_kernel=NULL   # not used below: frees the default kernel (N x N) during the gene loop

  # GxE multipliers: the numeric condition (one column), or indicators of the non-reference levels;
  # genotype x column j is the interaction term j (as in model.matrix(~ geno:condition))
  # A single condition (e.g. a constant) gives a plain eQTL map: no interaction terms.
  one_condition = length(unique(anno$condition)) == 1
  C = if (one_condition) matrix(0, nrow(anno), 0) else
    if (is.numeric(anno$condition)) cbind(anno$condition) else model.matrix(~ condition, data=anno)[, -1, drop=F]
  df_interact = ncol(C)

  no_geno = if (one_condition) matrix(0, nrow(anno), 0) else model.matrix( ~ condition, data=anno)
  if (!is.null(covariates)) {
    if (is.data.frame(covariates)) covariates = model.matrix(~ ., data = covariates)
    covariates = as.matrix(covariates)
    if (nrow(covariates) != nrow(anno) || !is.numeric(covariates) || anyNA(covariates))
      stop("covariates must be numeric, without NA, with one row per sample (as anno).", call. = FALSE)
    no_geno = cbind(no_geno, covariates)
  }
  no_geno = no_geno[, apply(no_geno, 2, sd) > 0, drop=F] # intercept: removed by the contrasts
  X0 = standardize_cols(no_geno)
  UX0 = crossprod(U, X0)

  # condition varying within individual => GxE columns can fit the kernel's null directions
  within = any(tapply(as.character(anno$condition), anno$individual, function(x) length(unique(x)) > 1))
  k_anc = if (la_mode) length(local_ancestry$dosage) else 1
  check_identifiable(eig, X0, need = 1 + if (within) (2 * k_anc - 1) * df_interact else 0)

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
  res = foreach(gene=genes, seed=seeds, .errorhandling=errorhandling) %dopar% {

    set.seed(seed)
    if (!is.null(checkpoint_dir)){
      check_fn = file.path(checkpoint_dir, paste0(gene, ".txt.gz"))
      if (file.exists(check_fn)) {
        return(read_tsv(check_fn))
      }
    }
    
    if (la_mode) {
      rows = cis_rows(geneloc[gene,"chr"], geneloc[gene,"left"]-cisdist, geneloc[gene,"right"]+cisdist)
      if (!length(rows)) return(NULL)
      D = lapply(local_ancestry$dosage, function(m) m[rows,,drop=F])
      H = lapply(local_ancestry$hapcount, function(m) m[rows,,drop=F])
      if (permutation_approach=="permute") {   # same shuffle of individuals for every ancestry
        perm = sample(ncol(D[[1]]))
        D = lapply(D, function(m) { colnames(m) = colnames(m)[perm]; m }); H = lapply(H, function(m) { colnames(m) = colnames(m)[perm]; m })
      }
      Uy = drop(crossprod(U, as.numeric(input[gene,])))
      gene_results = la_gene(gene, D, H, anno, C, U, lambda, Uy, UX0, permutation_approach=="boot", maf_threshold)
      if (!is.null(checkpoint_dir) && !is.null(gene_results)) gene_results %>% write_tsv(check_fn)
      return(gene_results)
    }

    geno_here = if (is.matrix(genotype)) {
      rows = cis_rows(geneloc[gene,"chr"], geneloc[gene,"left"]-cisdist, geneloc[gene,"right"]+cisdist)
      if (length(rows)) genotype[rows,,drop=F] else NULL
    } else get_cis_geno(genotype, 
                              chrom = geneloc[gene,"chr"], 
                              left = geneloc[gene,"left"]-cisdist, 
                              right = geneloc[gene,"right"]+cisdist, 
                              snploc = snploc,
                              genome_build = genome_build)
    if (is.null(geno_here)) return(NULL)
    
    imp_geno = if (anyNA(geno_here)) easy_impute(geno_here) else {   # nothing to impute: rounded, as easy_impute does
      g = round(geno_here); storage.mode(g) = "integer"; g }
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
    keep = rowSums(imp_geno) >= 5 & apply(G, 2, sd) > 0            # >= 5 ALT alleles, varies across samples
    if (!any(keep)) return(NULL)
    G = G[, keep, drop=F]
    Uy = drop(crossprod(U, as.numeric(input[gene,])))

    if (debug) cat("Fitting", ncol(G), "SNPs...\n")
    fit = .Call(suez_fit_gene, U, lambda, Uy, UX0, G, C, permutation_approach=="boot")   # src/lmm.cpp
    m = fit$res
    gene_results = data.frame(gene=gene, cis_snp=colnames(G), l0=fit$l0, l_geno=m[,1], l_interact=m[,2],
                              stringsAsFactors=F)
    gene_results$p_geno = pchisq(pmax(0, 2*(gene_results$l_geno-gene_results$l0)), 1, lower.tail=F)
    gene_results$df_interact = as.integer(m[,3])
    gene_results$p_interact = pchisq(pmax(0, 2*(gene_results$l_interact-gene_results$l_geno)), m[,3], lower.tail=F)
    gene_results$beta_geno = m[,4]
    gene_results$se_geno = m[,5]
    if (permutation_approach=="boot") {
      # parametric bootstrap under the genotype (no-interaction) model, one draw per SNP
      gene_results$l_boot_geno = m[,6]
      gene_results$l_boot_interact = m[,7]
    }

    if (!is.null(checkpoint_dir)){
      cat("Saving results to ", checkpoint_dir, "\n")
      gene_results %>% write_tsv(check_fn)
    }

    gene_results
  }
  bind_rows(res)
}
