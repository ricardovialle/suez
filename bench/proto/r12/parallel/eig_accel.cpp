// Symmetric eigen-decomposition: Eigen's solver (as suez_eigen_sym) vs Apple Accelerate's LAPACK dsyevd / dsyevr
// (multithreaded on Apple silicon; VECLIB_MAXIMUM_THREADS caps it). macOS only:
//   Sys.setenv(PKG_LIBS = "-framework Accelerate"); Rcpp::sourceCpp("eig_accel.cpp")
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using Eigen::MatrixXd; using Eigen::VectorXd;
extern "C" {
void dsyevd_(const char* jobz, const char* uplo, const int* n, double* a, const int* lda, double* w, double* work,
             const int* lwork, int* iwork, const int* liwork, int* info);
void dsyevr_(const char* jobz, const char* range, const char* uplo, const int* n, double* a, const int* lda,
             const double* vl, const double* vu, const int* il, const int* iu, const double* abstol, int* m,
             double* w, double* z, const int* ldz, int* isuppz, double* work, const int* lwork, int* iwork,
             const int* liwork, int* info);
}

// [[Rcpp::export]]
Rcpp::List eig_eigen(const Eigen::Map<MatrixXd> M) {
  Eigen::SelfAdjointEigenSolver<MatrixXd> es(M);
  return Rcpp::List::create(Rcpp::Named("values") = VectorXd(es.eigenvalues().reverse()),
                            Rcpp::Named("vectors") = MatrixXd(es.eigenvectors().rowwise().reverse()));
}

// [[Rcpp::export]]
Rcpp::List eig_accel(Rcpp::NumericMatrix M, bool rrr) {
  const int n = M.nrow(); Rcpp::NumericMatrix A = Rcpp::clone(M); Rcpp::NumericVector w(n); int info = 0;
  if (!rrr) {
    int lwork = -1, liwork = -1, iw; double wk;
    dsyevd_("V", "L", &n, A.begin(), &n, w.begin(), &wk, &lwork, &iw, &liwork, &info);
    lwork = (int)wk; liwork = iw; std::vector<double> work(lwork); std::vector<int> iwork(liwork);
    dsyevd_("V", "L", &n, A.begin(), &n, w.begin(), work.data(), &lwork, iwork.data(), &liwork, &info);
    if (info) Rcpp::stop("dsyevd info %d", info);
    MatrixXd V = Eigen::Map<MatrixXd>(A.begin(), n, n);
    return Rcpp::List::create(Rcpp::Named("values") = VectorXd(Eigen::Map<VectorXd>(w.begin(), n).reverse()),
                              Rcpp::Named("vectors") = MatrixXd(V.rowwise().reverse()));
  }
  Rcpp::NumericMatrix Z(n, n); std::vector<int> isuppz(2 * n); int m, lwork = -1, liwork = -1, iw; double wk, z0 = 0, tol = 0; int i0 = 0;
  dsyevr_("V", "A", "L", &n, A.begin(), &n, &z0, &z0, &i0, &i0, &tol, &m, w.begin(), Z.begin(), &n, isuppz.data(), &wk, &lwork, &iw, &liwork, &info);
  lwork = (int)wk; liwork = iw; std::vector<double> work(lwork); std::vector<int> iwork(liwork);
  dsyevr_("V", "A", "L", &n, A.begin(), &n, &z0, &z0, &i0, &i0, &tol, &m, w.begin(), Z.begin(), &n, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
  if (info) Rcpp::stop("dsyevr info %d", info);
  MatrixXd V = Eigen::Map<MatrixXd>(Z.begin(), n, n);
  return Rcpp::List::create(Rcpp::Named("values") = VectorXd(Eigen::Map<VectorXd>(w.begin(), n).reverse()),
                            Rcpp::Named("vectors") = MatrixXd(V.rowwise().reverse()));
}
