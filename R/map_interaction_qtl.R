

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
  apply(input, 1, qqnorm_no_plot) %>% t()
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
#' @param genotype [SNPs x individuals] matrix of genotypes (can have some NAs) OR a Bioconductor `TabixFile`` object
#' @param geneloc [genes x 4] data.frame with columns (geneid, chr, left, right). The intersect of geneloc$geneid and rownames(input) will be tested.
#' @param snploc [SNPs x 3] data.frame with columns (snpid, chr, pos) where snpid's map to rows of `genotype`
#' @param anno [samples x 2] data.frame with columns (individual, condition) where individual must correspond to column names of `genotype`
#' @param sample_kernel [samples x samples] covariance matrix learned in suez step 1. If you don't want to control for latent confounders leave this as NULL to use the matrix representing which samples are from the same individual.
#' @param normalization_approach One of (qq,log,linear) determining how to normalize within each gene. We recommend qq (quantile normalization to normal) to ensure the assumptions of the test hold.
#' @param permutation_approach One of ("none","permute","boot"). We recommend boot (parametric bootstrap) since the permutation test is not really valid for interaction effect testing.
#' @param cisdist How far from the gene boundaries to look for SNPs to test
#' @param checkpoint_dir An optional directory to store per gene results. If mapping crashes these results will be reused to save time.
#' @param debug Controls whether you see progress messages and whether errors are suppressed.
#' @param genome_build Only required if `genotype` is `TabixFile`. Default is "GRCh38"
#' @param maf_threshold Minimum MAF to analyze variants. 
#'
#' @details Model: y = X b + e, e ~ N(0, sigma2_k * K + sigma2 * I), K = `sample_kernel`, fitted by
#' maximum likelihood on the N-1 contrasts orthogonal to the intercept (REML for the intercept only),
#' so a kernel with K 1 = 0 (a GRM from centred genotypes) is handled. Each fit maximises the
#' likelihood exactly: b and the total variance have closed forms given h = sigma2_k / (sigma2_k +
#' sigma2), and h is found by grid and golden-section search on [0, 1]. Failed fits give NA. Columns are standardized
#' before fitting (the likelihood is invariant). The function stops if the kernel is not positive
#' semi-definite, or if the condition design leaves no kernel null direction for sigma2 (e.g. the
#' default kernel with condition fitting every within-individual contrast).
#'
#' @return data.frame, one row per gene x SNP: `l0`, `l_geno`, `l_interact` (log-likelihoods of the
#' contrasts under the no-genotype, genotype and genotype x condition models), `p_geno` (LRT, df 1),
#' `df_interact` and `p_interact` (LRT of `l_interact` vs `l_geno`), `beta_geno` and `se_geno` (GLS
#' genotype effect per ALT allele and its SE at the fitted variances), and for
#' `permutation_approach = "boot"` `l_boot_geno` and `l_boot_interact` (one parametric bootstrap draw
#' under the genotype model: 2 * (l_boot_interact - l_boot_geno) is a null draw of the interaction LR).
#' @importFrom rstan optimizing
#' @import doMC
#' @useDynLib suez, .registration = TRUE
#' @export
map_interaction_qtl = function(input, genotype, geneloc, snploc, anno, sample_kernel=NULL, 
                               normalization_approach="qq", permutation_approach="boot", cisdist=1e5, 
                               checkpoint_dir=NULL, debug=F, genome_build= "GRCh38", maf_threshold = 0) {

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
  rownames(geneloc)=geneloc$geneid

  if (is.null(sample_kernel)) {
    sample_kernel=outer(anno$individual,anno$individual,"==")
    class(sample_kernel)="numeric"
  }

  cat("Performing initial eigendecomposition...\n")
  eig=kernel_contrast_eigen(sample_kernel) # N-1 contrasts orthogonal to the intercept
  U=eig$vectors
  lambda=eig$values

  # GxE design for a genotype vector (numeric condition: one column; factor: one per non-reference level)
  interact_design=function(geno) {
    interact=model.matrix(~geno:condition,data=cbind(anno, geno=geno))
    if (is.numeric(anno$condition)) interact[,2,drop=F] else interact[,3:ncol(interact),drop=F]
  }
  df_interact=ncol(interact_design(rep(1, nrow(anno))))

  no_geno = model.matrix( ~ condition, data=anno)
  no_geno = no_geno[, apply(no_geno, 2, sd) > 0, drop=F] # intercept: removed by the contrasts
  X0 = standardize_cols(no_geno)
  UX0 = crossprod(U, X0)

  # condition varying within individual => GxE columns can fit the kernel's null directions
  within = any(tapply(as.character(anno$condition), anno$individual, function(x) length(unique(x)) > 1))
  check_identifiable(eig, X0, need = 1 + if (within) df_interact else 0)

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

  foreach(gene=genes, seed=seeds, .errorhandling=errorhandling, .combine = bind_rows) %dopar% {

    set.seed(seed)
    if (!is.null(checkpoint_dir)){
      check_fn = file.path(checkpoint_dir, paste0(gene, ".txt.gz"))
      if (file.exists(check_fn)) {
        return(read_tsv(check_fn))
      }
    }
    
    geno_here = get_cis_geno(genotype, 
                              chrom = geneloc[gene,"chr"], 
                              left = geneloc[gene,"left"]-cisdist, 
                              right = geneloc[gene,"right"]+cisdist, 
                              snploc = snploc,
                              genome_build = genome_build)
    if (is.null(geno_here)) return(NULL)
    
    imp_geno = easy_impute(geno_here)
    if (maf_threshold > 0) {
      mafs = rowMeans(geno_here, na.rm=T)/2 # calculate on unimputed
      imp_geno = imp_geno[mafs > maf_threshold,,drop=F]
      if (nrow(imp_geno)==0) return(NULL)
    }
    
    cis_snps = rownames(imp_geno)
    cat(gene,length(cis_snps)," cis snps\n")
    
    if (permutation_approach=="permute") colnames(imp_geno)=colnames(imp_geno)[ sample(ncol(imp_geno),ncol(imp_geno)) ]

    y=input[gene,] %>% as.numeric
    Uy=crossprod(U, y) %>% as.numeric

    if (debug) cat("Fitting no genotype model...\n")
    fit_no_geno=fit_lmm(Uy, UX0, lambda)
    if (is.null(fit_no_geno)) stop("no-genotype model did not converge for ", gene)

    gene_results = foreach(cis_snp=cis_snps, .errorhandling=errorhandling, .combine = bind_rows) %do% {

      geno=imp_geno[cis_snp,anno$individual]
      if (sum(imp_geno[cis_snp,]) < 5.0 || sd(geno) == 0) { return(NULL) }

      UX1=cbind(UX0, crossprod(U, standardize_cols(cbind(geno))))
      interact=interact_design(geno)
      interact=interact[, apply(interact, 2, sd) > 0, drop=F]
      UX2=cbind(UX1, crossprod(U, standardize_cols(interact)))

      lrt = function(Uy, f0) {
        fit_geno=fit_nested(Uy, UX1, lambda, f0)
        fit_interact=if (is.null(fit_geno) || ncol(interact)==0) NULL else fit_nested(Uy, UX2, lambda, fit_geno)
        list( fit_geno=fit_geno, fit_interact=fit_interact )
      }
      value=function(fit) if (is.null(fit)) NA_real_ else fit$value

      lrt_true=lrt( Uy, fit_no_geno )
      fit_geno=lrt_true$fit_geno
      fit_interact=lrt_true$fit_interact

      res=data.frame(gene=gene, cis_snp=cis_snp, l0=fit_no_geno$value, l_geno=value(fit_geno),
                     l_interact=value(fit_interact), stringsAsFactors=F  )
      res$p_geno=pchisq(pmax(0, 2*(res$l_geno-res$l0)), 1, lower.tail=F)
      res$df_interact=ncol(interact)
      res$p_interact=pchisq(pmax(0, 2*(res$l_interact-res$l_geno)), ncol(interact), lower.tail=F)
      bse=if (is.null(fit_geno)) c(NA, NA) else gls_coef(fit_geno, UX1, Uy, lambda) / sd(geno) # per allele
      res$beta_geno=bse[1]
      res$se_geno=bse[2]

      if (permutation_approach=="boot" && !is.null(fit_geno)) {
        # parametric bootstrap under the genotype (no-interaction) model, in the rotated space
        o = fit_geno$par$sigma2_k * lambda + fit_geno$par$sigma2
        Uy_boot = sqrt(o) * rnorm(length(o)) + drop(UX1 %*% fit_geno$par$beta)
        fit_boot_0 = fit_lmm(Uy_boot, UX0, lambda)
        lrt_boot = if (is.null(fit_boot_0)) list() else lrt( Uy_boot, fit_boot_0 )
        res$l_boot_geno=value(lrt_boot$fit_geno)
        res$l_boot_interact=value(lrt_boot$fit_interact)
      }

      res
    }
    if (!is.null(checkpoint_dir)){
      cat("Saving results to ", checkpoint_dir, "\n")
      gene_results %>% write_tsv(check_fn)
    }

    gene_results
  }
}

