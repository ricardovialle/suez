// Symmetric eigen-decomposition (ascending values, vectors) with Accelerate's LAPACK (dsyevd / dsyevr) or Eigen.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
extern "C" {
void dsyevd_(const char*, const char*, const int*, double*, const int*, double*, double*, const int*, int*, const int*, int*);
void dsyevr_(const char*, const char*, const char*, const int*, double*, const int*, const double*, const double*, const int*, const int*,
             const double*, int*, double*, double*, const int*, int*, double*, const int*, int*, const int*, int*);
}
// [[Rcpp::export]]
Rcpp::List eig_acc(Rcpp::NumericMatrix A0, int method) {   // 0 = Eigen, 1 = dsyevd, 2 = dsyevr
  const int n = A0.nrow(); Eigen::MatrixXd A = Rcpp::as<Eigen::MatrixXd>(A0); Eigen::VectorXd w(n); int info = 0;
  if (method == 0) { Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A); return Rcpp::List::create(Rcpp::Named("values") = es.eigenvalues(), Rcpp::Named("vectors") = es.eigenvectors()); }
  if (method == 1) {
    int lwork = -1, liwork = -1, iw; double wq;
    dsyevd_("V", "L", &n, A.data(), &n, w.data(), &wq, &lwork, &iw, &liwork, &info);
    lwork = (int)wq; liwork = iw; std::vector<double> work(lwork); std::vector<int> iwork(liwork);
    dsyevd_("V", "L", &n, A.data(), &n, w.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
    return Rcpp::List::create(Rcpp::Named("values") = w, Rcpp::Named("vectors") = A, Rcpp::Named("info") = info);
  }
  Eigen::MatrixXd Z(n, n); std::vector<int> isuppz(2 * n); int m, lwork = -1, liwork = -1, iw; double wq, vl = 0, vu = 0, abstol = 0; int il = 0, iu = 0;
  dsyevr_("V", "A", "L", &n, A.data(), &n, &vl, &vu, &il, &iu, &abstol, &m, w.data(), Z.data(), &n, isuppz.data(), &wq, &lwork, &iw, &liwork, &info);
  lwork = (int)wq; liwork = iw; std::vector<double> work(lwork); std::vector<int> iwork(liwork);
  dsyevr_("V", "A", "L", &n, A.data(), &n, &vl, &vu, &il, &iu, &abstol, &m, w.data(), Z.data(), &n, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
  return Rcpp::List::create(Rcpp::Named("values") = w, Rcpp::Named("vectors") = Z, Rcpp::Named("info") = info);
}
