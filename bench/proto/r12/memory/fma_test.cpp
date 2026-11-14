#include <Rcpp.h>
#include <cmath>
#pragma clang fp contract(off)
// [[Rcpp::export]]
Rcpp::NumericVector tdot(Rcpp::NumericMatrix A, Rcpp::NumericVector x, bool fma) {
  int n = A.nrow(), m = A.ncol(); Rcpp::NumericVector z(m);
  for (int j = 0; j < m; ++j) { double t = 0; const double* a = &A(0, j);
    for (int i = 0; i < n; ++i) t = fma ? std::fma(a[i], x[i], t) : t + a[i] * x[i]; z[j] = t; }
  return z;
}
// [[Rcpp::export]]
Rcpp::NumericVector ndot(Rcpp::NumericMatrix A, Rcpp::NumericVector x, bool fma) {
  int n = A.nrow(), m = A.ncol(); Rcpp::NumericVector y(n);
  for (int j = 0; j < m; ++j) { double t = x[j]; const double* a = &A(0, j);
    for (int i = 0; i < n; ++i) y[i] = fma ? std::fma(t, a[i], y[i]) : y[i] + t * a[i]; }
  return y;
}
