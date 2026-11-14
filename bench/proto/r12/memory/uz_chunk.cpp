// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using Eigen::MatrixXd;
// base rotate's repeated-row branch vs the same with UZ built and multiplied in row blocks of the result
// [[Rcpp::export]]
Rcpp::List uz(const Eigen::Map<MatrixXd> U, const Eigen::Map<MatrixXd> M, int rb) {
  const int N = M.rows(), S = M.cols();
  const MatrixXd Mt = M.transpose();
  std::vector<int> ord(N); for (int i = 0; i < N; ++i) ord[i] = i;
  auto less = [&](int a, int b) { const double *x = Mt.col(a).data(), *y = Mt.col(b).data();
    for (int j = 0; j < S; ++j) if (x[j] != y[j]) return x[j] < y[j]; return false; };
  std::sort(ord.begin(), ord.end(), less);
  std::vector<int> cls(N); int nc = 0;
  for (int i = 0; i < N; ++i) cls[ord[i]] = (i > 0 && !less(ord[i - 1], ord[i])) ? nc - 1 : nc++;
  MatrixXd UZ = MatrixXd::Zero(U.cols(), nc), Mu(nc, S);
  for (int i = 0; i < N; ++i) { UZ.col(cls[i]) += U.row(i).transpose(); Mu.row(cls[i]) = M.row(i); }
  MatrixXd R1 = UZ * Mu; UZ.resize(0, 0);
  MatrixXd R2(U.cols(), S);
  for (int r0 = 0; r0 < U.cols(); r0 += rb) {
    const int h = std::min(rb, (int)U.cols() - r0);
    MatrixXd Z = MatrixXd::Zero(h, nc);
    for (int i = 0; i < N; ++i) Z.col(cls[i]) += U.row(i).segment(r0, h).transpose();
    R2.middleRows(r0, h) = Z * Mu;
  }
  return Rcpp::List::create(R1, R2, nc);
}
