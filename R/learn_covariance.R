#' Function to learn covariance structure across samples using latent factors
#'
#' @param input [genes x samples] Normalized (e.g. quantile normalized) matrix of log2 gene expression
#' @param anno [samples x 2] data.frame with `individual` and `condition` columns. `condition` may be a factor/character (one effect per level) or numeric (one slope). If a kinship matrix is provided the row and column names must correspond to the values used in the `individual` column.
#' @param P Number of latent factors to model. suez can prune out unnecessary factors (in principle).
#' @param kinship_matrix Representing relatedness between individuals. Optional.
#' @param same_ind_initial_weight Initialization.
#' @param kinship_matrix_initial_weight Initialization.
#' @param ... passed to rstan::optimizing
#'
#' @details Model, per gene: y = X b + e, e ~ N(0, s_noise I + s_ind SameInd + s_kin Kinship +
#' sum_p s_p x_p x_p'), X = model.matrix(~condition). The variances and latent factors (unit vectors
#' x_p) are fitted by maximum likelihood of the N - rank(X) contrasts orthogonal to X (REML): the
#' condition effects are removed exactly, without treating the regression residuals (which have no
#' variance along X) as N-dimensional data. A same-condition kernel lies in the span of X, so it has no
#' variance left to model and is not a component.
#'
#' @return list: `cov_suez` [samples x samples], the latent-factor, same-individual and kinship parts
#' of the covariance (no isotropic noise: pass it as `sample_kernel` to [map_interaction_qtl]), and
#' `optimized_fit`, the rstan::optimizing result (factors `x_adjustable` are in the contrast basis).
#' @importFrom stats model.matrix
#' @import doMC
#' @export
learn_covariance = function(input, anno, kinship_matrix=NULL, P=20, same_ind_initial_weight=0.01, kinship_matrix_initial_weight=0.01, ... )
{

  self_outer=function(g) outer(g,g)

  condition = if (is.numeric(anno$condition)) anno$condition else as.factor(anno$condition)
  X = model.matrix( ~ condition)
  q = qr(X)
  # orthonormal basis of the contrasts orthogonal to the condition design
  Q = qr.Q(q, complete = TRUE)[, -seq_len(q$rank), drop = FALSE]
  resi = input %*% Q   # [genes x (N - rank X)]
  rotate = function(K) crossprod(Q, K %*% Q)

  svd_ge=svd(resi)

  # Latent factors initialized using SVD
  x_adjustable=foreach(p=1:P) %do% { svd_ge$v[,p] }
  sinit=svd_ge$d^2 / nrow(input)

  same_ind=outer(anno$individual, anno$individual, "==") * 1

  variance_components=list( isotropic_noise=diag(ncol(resi)),
                           same_ind=rotate(same_ind))
  # residual after removing first P PCs (used to estimate noise variance)
  rr = resi - svd_ge$u[,1:P] %*% diag(svd_ge$d[1:P], P) %*% t(svd_ge$v[,1:P])
  variance_component_weights=c( mean(rr^2), same_ind_initial_weight )

  if (!is.null(kinship_matrix)) {
    kin = kinship_matrix[ anno$individual, anno$individual ]
    variance_components$kinship_matrix=rotate(kin)
    variance_component_weights=c(variance_component_weights, kinship_matrix_initial_weight)
  }

  dat=list(YtY=crossprod(resi),
           N=ncol(resi),
           G=nrow(resi),
           P=length(variance_components),
           x=variance_components,
           P_adjustable=length(x_adjustable))

  init=list( s_adjustable=sinit[1:P], s=variance_component_weights, x_adjustable=x_adjustable )

  optimized_fit=rstan::optimizing(stanmodels$suez_step_1,
                           data=dat,
                           init=init,
                           as_vector=F,
                           ...)
  if (optimized_fit$return_code != 0)
    warning("rstan::optimizing did not report convergence (return code ", optimized_fit$return_code,
            "); consider more iterations (iter = ...).", call. = FALSE)

  names(optimized_fit$par$s) = names(variance_components)

  # factors back to sample space (Q x is still a unit vector)
  cov_suez = Reduce("+", foreach(p=1:dat$P_adjustable) %do% {
    optimized_fit$par$s_adjustable[p] * self_outer(drop(Q %*% optimized_fit$par$x_adjustable[p,]))
    } )

  cov_suez = cov_suez + optimized_fit$par$s["same_ind"] * same_ind
  if (!is.null(kinship_matrix)) {
    cov_suez = cov_suez + optimized_fit$par$s["kinship_matrix"] * kin
  }

  list( cov_suez = cov_suez, optimized_fit = optimized_fit )
}
