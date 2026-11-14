// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <chrono>
using Eigen::MatrixXd; using Eigen::VectorXd;
static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
// [[Rcpp::export]]
Rcpp::NumericVector mw_bench(int R, int P2, int iters) {
  MatrixXd MM = MatrixXd::Random(R, P2), AB = MatrixXd::Random(R, 3), MW(P2, 3), MMt = MM.transpose();
  double t[6] = {0}; double chk = 0;
  for (int v = 0; v < 6; ++v) {
    double t0 = now();
    for (int it = 0; it < iters; ++it) {
      AB(0, 0) += 1e-9;
      if (v == 0) MW.noalias() = MM.transpose() * AB;
      else if (v == 1) { MW.col(0).noalias() = MM.transpose() * AB.col(0); MW.middleCols(1, 2).noalias() = MM.transpose() * AB.middleCols(1, 2); }
      else if (v == 2) MW.noalias() = MM.transpose().lazyProduct(AB);
      else if (v == 3) { for (int q = 0; q < P2; ++q) for (int j = 0; j < 3; ++j) MW(q, j) = MM.col(q).dot(AB.col(j)); }
      else if (v == 4) MW.noalias() = MMt * AB;
      else if (v == 5) MW.col(0).noalias() = MM.transpose() * AB.col(0);
      chk += MW(0, 0);
    }
    t[v] = (now() - t0) / iters * 1e6;
  }
  Rcpp::NumericVector out(t, t + 6); out.attr("names") = Rcpp::CharacterVector::create("gemm3", "gemv+gemm2", "lazy", "dots", "rowmajor_gemm", "gemv_only");
  return out;
}
