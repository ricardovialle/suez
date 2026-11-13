data {
  int<lower=0> N; // individuals
  int<lower=0> G; // genes
  int<lower=0> P; // covariates
  int<lower=0> P_adjustable; // covariates
  array[P] matrix[N,N] x;
  matrix[N,N] YtY; // sum over genes of y_g y_g' (sufficient statistic)
}
transformed data {
  // fixed components as columns, so that their weighted sum is one matrix-vector product: the autodiff
  // tape then holds O(N^2) entries per gradient instead of O(N^2) per component (memory, see bench)
  matrix[N * N, P] xflat;
  for (p in 1:P) xflat[:, p] = to_vector(x[p]);
}
parameters {
  array[P] real<lower=0> s;
  array[P_adjustable] real<lower=0> s_adjustable;
  array[P_adjustable] unit_vector[N] x_adjustable;
}
model {
  matrix[N, P_adjustable] xa;
  matrix[N,N] sigma;
  matrix[N,N] L;
  for (p in 1:P_adjustable)
    xa[:, p] = sqrt(s_adjustable[p]) * x_adjustable[p];
  // = sum_p s[p] x[p] + sum_p s_adjustable[p] x_adjustable[p] x_adjustable[p]'
  sigma = to_matrix(xflat * to_vector(s), N, N) + tcrossprod(xa);
  L = cholesky_decompose(sigma);
  // = sum_g multi_normal_cholesky_lupdf(y[g] | 0, L), without a pass over genes
  target += - G * sum(log(diagonal(L))) - 0.5 * sum(YtY .* crossprod(mdivide_left_tri_low(L, diag_matrix(rep_vector(1.0, N)))));
}
