// Exact ML fit of Uy = UX b + e, e ~ N(0, diag(sigma2_k * lambda + sigma2)), for map_interaction_qtl.
// Writing sigma2_k = s h, sigma2 = s (1 - h), b and s have closed forms for each h in [0, 1], so the
// fit is a 1-D maximisation over h (grid, then golden section), boundaries included.
#include <RcppEigen.h>
#include <cmath>
#include <limits>

using Eigen::ArrayXd; using Eigen::MatrixXd; using Eigen::VectorXd;

struct ProfileLmm {
  const Eigen::Map<MatrixXd> X; const Eigen::Map<VectorXd> y; const Eigen::Map<VectorXd> lam;
  VectorXd b; double s = NAN;

  // profile log-likelihood at h (no 2 pi constant); sets b and s
  double operator()(double h) {
    const int n = y.size();
    ArrayXd d = h * lam.array() + (1 - h);
    if ((d <= 0).any()) return -std::numeric_limits<double>::infinity();
    ArrayXd w = d.inverse();
    VectorXd r = y;
    if (X.cols()) {
      MatrixXd Xw = X.array().colwise() * w;
      b = (Xw.transpose() * X).ldlt().solve(Xw.transpose() * y);
      r -= X * b;
    } else b.resize(0);
    s = (r.array().square() * w).sum() / n;
    if (!(s > 0)) return std::numeric_limits<double>::infinity();   // exact fit: unbounded
    return -0.5 * n - 0.5 * (n * std::log(s) + d.log().sum());
  }
};

extern "C" SEXP suez_profile_lmm(SEXP UX, SEXP Uy, SEXP lambda) {
  BEGIN_RCPP
  ProfileLmm f{Rcpp::as<Eigen::Map<MatrixXd>>(UX), Rcpp::as<Eigen::Map<VectorXd>>(Uy),
               Rcpp::as<Eigen::Map<VectorXd>>(lambda)};
  const int ng = 51;
  std::vector<double> v(ng); int best = 0;
  for (int i = 0; i < ng; ++i) { v[i] = f(i / (ng - 1.0)); if (v[i] > v[best]) best = i; }
  double hbest = best / (ng - 1.0), vbest = v[best];
  // golden section on the bracket around the best grid point
  double a = std::max(0, best - 1) / (ng - 1.0), c = std::min(ng - 1, best + 1) / (ng - 1.0);
  const double g = 0.5 * (std::sqrt(5.0) - 1);
  double x1 = c - g * (c - a), x2 = a + g * (c - a), f1 = f(x1), f2 = f(x2);
  while (c - a > 1e-10) {
    if (f1 > f2) { c = x2; x2 = x1; f2 = f1; x1 = c - g * (c - a); f1 = f(x1); }
    else         { a = x1; x1 = x2; f1 = f2; x2 = a + g * (c - a); f2 = f(x2); }
  }
  if (f1 > vbest) { hbest = x1; vbest = f1; }
  if (f2 > vbest) { hbest = x2; vbest = f2; }
  vbest = f(hbest);   // refresh b and s at the optimum
  return Rcpp::List::create(Rcpp::Named("value") = vbest, Rcpp::Named("sigma2_k") = f.s * hbest,
                            Rcpp::Named("sigma2") = f.s * (1 - hbest), Rcpp::Named("beta") = f.b);
  END_RCPP
}
