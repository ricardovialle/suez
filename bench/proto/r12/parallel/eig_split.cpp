// Where Eigen's SelfAdjointEigenSolver spends its time, to bound what Eigen's OpenMP (which parallelises only
// its GEMM-based products) could gain: tridiagonalisation (level-2 rank-2 updates: not parallel), forming Q
// (blocked Householder, GEMM: parallel with -fopenmp), implicit QR on the tridiagonal applying Givens rotations
// to Q (serial).
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <chrono>
using Eigen::MatrixXd; using Eigen::VectorXd;
static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
// [[Rcpp::export]]
Rcpp::NumericVector eig_split(const Eigen::Map<MatrixXd> A) {
  const int n = A.rows(); double t0 = now();
  MatrixXd M = A; VectorXd d(n), e(n - 1), hc(n - 1);
  Eigen::internal::tridiagonalization_inplace(M, d, e, hc, false);    // reduction only
  double t1 = now();
  MatrixXd M2 = A; Eigen::Tridiagonalization<MatrixXd> tr(M2); MatrixXd Q = tr.matrixQ();
  double t2 = now();
  Eigen::SelfAdjointEigenSolver<MatrixXd> es(A);
  double t3 = now();
  return Rcpp::NumericVector::create(Rcpp::Named("tridiag") = t1 - t0, Rcpp::Named("tridiag_plus_Q") = t2 - t1,
                                     Rcpp::Named("full") = t3 - t2);
}
