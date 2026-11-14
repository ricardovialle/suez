// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <chrono>
using Eigen::MatrixXd; using Eigen::VectorXd;
static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
// [[Rcpp::export]]
Rcpp::NumericVector gv(int R, int P2, int iters) {
  MatrixXd MM = MatrixXd::Random(R, P2), MMt = MM.transpose(); VectorXd a = VectorXd::Random(R), w1(P2), w2(P2);
  double t0 = now(), s = 0;
  for (int it = 0; it < iters; ++it) { a[0] += 1e-9; w1.noalias() = MM.transpose() * a; s += w1[0]; }
  double t1 = now();
  for (int it = 0; it < iters; ++it) { a[0] += 1e-9; w2.noalias() = MMt * a; s += w2[0]; }
  double t2 = now();
  return Rcpp::NumericVector::create((t1 - t0) / iters * 1e6, (t2 - t1) / iters * 1e6, (w1 - w2).norm() / w1.norm());
}
