// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
// [[Rcpp::export]]
Eigen::MatrixXd prodE(Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B) { return A.transpose() * B; }
