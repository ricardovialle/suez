#' The 'suez' package.
#' 
#' @description Maps expression QTLs and genotype x condition (response) eQTLs with linear mixed models that account for repeated samples, relatedness and hidden confounders, extending the PANAMA approach to repeat measurements and multiple conditions (with optional local-ancestry models for admixed samples)
#' 
#' @docType package
#' @name suez-package
#' @aliases suez
#' @useDynLib suez, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstantools rstan_config
#' @import doMC
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom magrittr %>%
#' @importFrom utils globalVariables
#' @importFrom stats lm.fit model.matrix optimize pchisq pt qnorm rnorm sd
#' 
#' @references 
#' * Fusi N, Stegle O, Lawrence ND. 2012. Joint modelling of confounding factors and prominent genetic regulators provides increased accuracy in genetical genomics studies. PLoS Computational Biology 8:e1002330.
#' 
#' * Knowles* DA, Burrows* CK, Blischak JD, Patterson KM, Serie DJ, Norton N, Ober C, Pritchard JK and Gilad Y (2018), "Determining the genetic basis of anthracycline-cardiotoxicity by molecular response QTL mapping in induced cardiomyocytes", eLife
#' DOI: https://doi.org/10.1371/journal.pcbi.1002330, PMID: 22241974
#' 
#' * Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.17.3. http://mc-stan.org
#' 
NULL

# foreach loop variables
utils::globalVariables(c("p", "gene", "cis_snp", "seed"))
