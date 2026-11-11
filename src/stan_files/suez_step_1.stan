data {
  int<lower=0> N; // individuals
  int<lower=0> G; // genes
  int<lower=0> P; // covariates
  int<lower=0> P_adjustable; // covariates
  array[P] matrix[N,N] x;
  matrix[N,N] YtY; // sum over genes of y_g y_g' (sufficient statistic)
}
parameters {
  array[P] real<lower=0> s;
  array[P_adjustable] real<lower=0> s_adjustable;
  array[P_adjustable] unit_vector[N] x_adjustable;
}
model {
  matrix[N,N] sigma;
  matrix[N,N] L;
  matrix[N,N] L_inv;
  sigma = rep_matrix(0.0, N, N);
  for (p in 1:P)
    sigma = sigma + s[p] * x[p];
  for (p in 1:P_adjustable)
    sigma = sigma + s_adjustable[p] * ( x_adjustable[p] * to_row_vector(x_adjustable[p]) );
  L=cholesky_decompose(sigma);
  L_inv = mdivide_left_tri_low(L, diag_matrix(rep_vector(1.0, N)));
  // = sum_g multi_normal_cholesky_lupdf(y[g] | 0, L), without a pass over genes
  target += - G * sum(log(diagonal(L))) - 0.5 * sum(YtY .* crossprod(L_inv));
}
