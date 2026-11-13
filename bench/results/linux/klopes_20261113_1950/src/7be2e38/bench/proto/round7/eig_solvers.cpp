// [[Rcpp::depends(RcppEigen)]]
// Round 7 prototype: symmetric eigen-decomposition with Eigen (RcppEigen, already a LinkingTo dependency)
// instead of R's eigen() (reference LAPACK dsyevr). Returns values in decreasing order, as eigen().
#include <RcppEigen.h>
// [[Rcpp::export]]
Rcpp::List eigen_sym(const Eigen::Map<Eigen::MatrixXd> A) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
  const int n = A.rows();
  Eigen::VectorXd v = es.eigenvalues().reverse();
  Eigen::MatrixXd V = es.eigenvectors().rowwise().reverse();
  return Rcpp::List::create(Rcpp::Named("values") = v, Rcpp::Named("vectors") = V);
}
