// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <chrono>
using Eigen::MatrixXd; using Eigen::VectorXd;
static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
static const double COLLINEAR = 1e-14;
struct PinvLDLT {
  Eigen::LDLT<MatrixXd> L; VectorXd d0;
  void compute(const MatrixXd& A) {
    L.compute(A); d0 = A.diagonal();
    const auto& T = L.transpositionsP(); for (int i = 0; i < T.size(); ++i) std::swap(d0[i], d0[T.coeff(i)]);
  }
  void solve(const VectorXd& rhs, VectorXd& x) const {   // x = A^+ rhs (no allocation once x has its size)
    x = rhs; const auto& T = L.transpositionsP(); const int m = x.size();
    for (int i = 0; i < m; ++i) std::swap(x[i], x[T.coeff(i)]);
    L.matrixL().solveInPlace(x);
    const VectorXd& D = L.vectorD();
    for (int i = 0; i < m; ++i) x[i] = std::fabs(D[i]) > COLLINEAR * std::fabs(d0[i]) && D[i] != 0 ? x[i] / D[i] : 0;
    L.matrixU().solveInPlace(x);
    for (int i = m - 1; i >= 0; --i) std::swap(x[i], x[T.coeff(i)]);
  }
  VectorXd solve(const VectorXd& rhs) const { VectorXd x; solve(rhs, x); return x; }
  std::vector<int> kept() const {   // original indices of the columns whose pivot is kept, increasing
    const auto& T = L.transpositionsP(); const VectorXd& D = L.vectorD(); const int m = D.size();
    std::vector<int> p(m), k; for (int i = 0; i < m; ++i) p[i] = i;
    for (int i = 0; i < m; ++i) std::swap(p[i], p[T.coeff(i)]);
    for (int i = 0; i < m; ++i) if (std::fabs(D[i]) > COLLINEAR * std::fabs(d0[i]) && D[i] != 0) k.push_back(p[i]);
    std::sort(k.begin(), k.end()); return k;
  }
};

// PinvLDLT's factorisation and pseudo-solve without Eigen's per-step GEMV overhead (small p, called per Newton
// step): Eigen's unblocked LDLT pivots on the largest *original* |diagonal| among the remaining columns (its
// diagonal is not updated before the search), so the same permutation is simulated first, then a right-looking
// LDLT (column axpys) on the permuted lower triangle. Same pivots, same keep/drop rule; rounding order differs.
struct FastLDLT {
  int m = 0; MatrixXd B; VectorXd D, d0; std::vector<int> perm; VectorXd x_;
  void compute(const MatrixXd& A) {
    m = A.rows(); perm.resize(m); D.resize(m); d0.resize(m); B.resize(m, m);
    for (int i = 0; i < m; ++i) { perm[i] = i; d0[i] = A(i, i); }
    for (int k = 0; k < m; ++k) {   // as Eigen: first index of the largest |d0| among positions k..m-1, then swap
      int bi = k; double bv = std::fabs(d0[k]);
      for (int i = k + 1; i < m; ++i) if (std::fabs(d0[i]) > bv) { bv = std::fabs(d0[i]); bi = i; }
      std::swap(perm[k], perm[bi]); std::swap(d0[k], d0[bi]);
    }
    for (int c = 0; c < m; ++c) for (int r = c; r < m; ++r) B(r, c) = A(perm[r], perm[c]);
    for (int k = 0; k < m; ++k) {
      const double dk = B(k, k); D[k] = dk; double* lk = &B(0, k);
      if (dk == 0) continue;   // Eigen: column left unscaled, and D = 0 removes it from every later update
      for (int i = k + 1; i < m; ++i) lk[i] /= dk;
      const double* __restrict lkr = lk;
      int j = k + 1;
      for (; j + 1 < m; j += 2) {   // two columns per pass: lk loaded once for both
        const double t0 = lkr[j] * dk, t1 = lkr[j + 1] * dk; double* __restrict b0 = &B(0, j); double* __restrict b1 = &B(0, j + 1);
        b0[j] -= t0 * lkr[j];
        for (int i = j + 1; i < m; ++i) { const double l = lkr[i]; b0[i] -= t0 * l; b1[i] -= t1 * l; }
      }
      for (; j < m; ++j) { const double t = lkr[j] * dk; double* __restrict bj = &B(0, j); for (int i = j; i < m; ++i) bj[i] -= t * lkr[i]; }
    }
  }
  void solve(const VectorXd& rhs, VectorXd& x) {   // x = A^+ rhs as PinvLDLT::solve
    x_.resize(m); for (int i = 0; i < m; ++i) x_[i] = rhs[perm[i]];
    for (int k = 0; k < m; ++k) { const double xk = x_[k]; const double* lk = &B(0, k); for (int i = k + 1; i < m; ++i) x_[i] -= lk[i] * xk; }
    for (int i = 0; i < m; ++i) x_[i] = std::fabs(D[i]) > COLLINEAR * std::fabs(d0[i]) && D[i] != 0 ? x_[i] / D[i] : 0;
    for (int k = m - 1; k >= 0; --k) { const double* lk = &B(0, k); double acc = x_[k]; for (int i = k + 1; i < m; ++i) acc -= lk[i] * x_[i]; x_[k] = acc; }
    x.resize(m); for (int i = 0; i < m; ++i) x[perm[i]] = x_[i];
  }
};


// [[Rcpp::export]]
Rcpp::NumericVector ldlt3(int p, int iters) {
  MatrixXd X = MatrixXd::Random(3 * p, p); MatrixXd A = X.transpose() * X; VectorXd b = VectorXd::Random(p), x, y;
  PinvLDLT P; FastLDLT F; double t0 = now(), s = 0;
  for (int it = 0; it < iters; ++it) { A(0, 0) += 1e-12; P.compute(A); s += P.L.vectorD()[0]; }
  double t1 = now();
  for (int it = 0; it < iters; ++it) { A(0, 0) += 1e-12; F.compute(A); s += F.D[0]; }
  double t2 = now();
  for (int it = 0; it < iters; ++it) { P.solve(b, x); s += x[0]; }
  double t3 = now();
  for (int it = 0; it < iters; ++it) { F.solve(b, y); s += y[0]; }
  double t4 = now(); P.solve(b, x);
  return Rcpp::NumericVector::create((t1 - t0) / iters * 1e6, (t2 - t1) / iters * 1e6, (t3 - t2) / iters * 1e6, (t4 - t3) / iters * 1e6, (x - y).norm() / x.norm());
}
