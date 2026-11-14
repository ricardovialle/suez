// Is Eigen's GEMM bitwise invariant to partitioning the output rows (or columns)? If so, a threaded U'G over
// row ranges of the result gives exactly the serial result at any thread count.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using Eigen::MatrixXd;
// [[Rcpp::export]]
Rcpp::NumericVector gemm_bits(const Eigen::Map<MatrixXd> A, const Eigen::Map<MatrixXd> B, int rchunk, int cchunk) {
  const MatrixXd full = A.transpose() * B;           // (A'B), as U'M in rotate()
  MatrixXd part(full.rows(), full.cols());
  for (int r0 = 0; r0 < full.rows(); r0 += rchunk) {
    const int h = std::min<int>(rchunk, full.rows() - r0);
    for (int c0 = 0; c0 < full.cols(); c0 += cchunk) {
      const int w = std::min<int>(cchunk, full.cols() - c0);
      part.block(r0, c0, h, w).noalias() = A.middleCols(r0, h).transpose() * B.middleCols(c0, w);
    }
  }
  // the UZ * Mu form (no transpose)
  const MatrixXd At = A.transpose(); const MatrixXd full2 = At * B; MatrixXd part2(full2.rows(), full2.cols());
  for (int r0 = 0; r0 < full2.rows(); r0 += rchunk) {
    const int h = std::min<int>(rchunk, full2.rows() - r0);
    part2.middleRows(r0, h).noalias() = At.middleRows(r0, h) * B;
  }
  return Rcpp::NumericVector::create((full - part).cwiseAbs().maxCoeff(), (full2 - part2).cwiseAbs().maxCoeff(),
                                     (full - full2).cwiseAbs().maxCoeff());
}
