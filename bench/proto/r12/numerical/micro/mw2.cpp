// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <chrono>
using Eigen::MatrixXd; using Eigen::VectorXd;
static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
static void axpy3(const MatrixXd& MMt, const MatrixXd& AB, MatrixXd& MW, int nc) {
  const int P2 = MMt.rows(), R = MMt.cols(); MW.setZero();
  double *w0 = MW.col(0).data(), *w1 = MW.col(1).data(), *w2 = MW.col(2).data();
  for (int r = 0; r < R; ++r) {
    const double* m = MMt.col(r).data(); const double a = AB(r, 0);
    if (nc == 1) { for (int q = 0; q < P2; ++q) w0[q] += a * m[q]; continue; }
    const double b = AB(r, 1), e = AB(r, 2);
    for (int q = 0; q < P2; ++q) { const double x = m[q]; w0[q] += a * x; w1[q] += b * x; w2[q] += e * x; }
  }
}
// [[Rcpp::export]]
Rcpp::NumericVector mw_bench2(int R, int P2, int iters) {
  MatrixXd MM = MatrixXd::Random(R, P2), AB = MatrixXd::Random(R, 3), MW(P2, 3), MW2(P2, 3), MMt = MM.transpose();
  double t[4] = {0}; double chk = 0;
  for (int v = 0; v < 4; ++v) {
    double t0 = now();
    for (int it = 0; it < iters; ++it) {
      AB(0, 0) += 1e-9;
      if (v == 0) MW.noalias() = MM.transpose() * AB;
      else if (v == 1) axpy3(MMt, AB, MW2, 3);
      else if (v == 2) MW.col(0).noalias() = MM.transpose() * AB.col(0);
      else axpy3(MMt, AB, MW2, 1);
      chk += MW(0, 0) + MW2(0, 0);
    }
    t[v] = (now() - t0) / iters * 1e6;
  }
  MW.noalias() = MM.transpose() * AB; axpy3(MMt, AB, MW2, 3);
  Rcpp::NumericVector out = {t[0], t[1], t[2], t[3], (MW - MW2).cwiseAbs().maxCoeff() / MW.cwiseAbs().maxCoeff()};
  out.attr("names") = Rcpp::CharacterVector::create("gemm3", "axpy3", "gemv1", "axpy1", "reldiff");
  return out;
}
// [[Rcpp::export]]
Rcpp::NumericVector ldlt_bench(int p, int iters) {
  MatrixXd X = MatrixXd::Random(3 * p, p); MatrixXd A = X.transpose() * X; VectorXd b = VectorXd::Random(p), x;
  Eigen::LDLT<MatrixXd> L; Eigen::LLT<MatrixXd> C; double t0 = now(), s = 0;
  for (int it = 0; it < iters; ++it) { A(0, 0) += 1e-12; L.compute(A); x = L.solve(b); s += x[0]; }
  double t1 = now();
  for (int it = 0; it < iters; ++it) { A(0, 0) += 1e-12; C.compute(A); x = C.solve(b); s += x[0]; }
  double t2 = now();
  return Rcpp::NumericVector::create((t1 - t0) / iters * 1e6, (t2 - t1) / iters * 1e6);
}
