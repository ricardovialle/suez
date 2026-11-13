// [[Rcpp::depends(RcppEigen)]]
// Prototype p2 (algorithm specialist): exact compression of the n rotated contrasts.
//
// Every quantity of the profile likelihood (grid and Newton) is a sum over contrasts i of f(lambda_i, h) times
// a quadratic form in the row z_i = [X_i, y_i]. Contrasts are grouped by eigenvalue:
//  * non-null groups: |lambda_i - c_g| <= 1e-12 c_g (c_g = group minimum), weight taken at c_g (relative error
//    <= 1e-12 in every weight, for every h);
//  * the null group (lambda <= 1e-12 max): eigen-decomposition noise around 0 (values 0 .. ~1e-14). Its
//    weights are expanded in delta_i = lambda_i - c_g: f(lambda_i) = sum_k f^(k)(c_g)/k! delta_i^k, k <= KN.
//    The truncation is checked at every h (rho = h max(delta) / d_g, rho^(KN+1) <= 1e-13), otherwise the
//    evaluation falls back to the exact per-contrast computation (always at h = 1 when lambda has zeros).
// Within a group (and Taylor order k) the rows are replaced by the R factor of the QR of the group's rows
// (scaled by delta^(k/2)): sum_i z_i z_i' = R'R and every residual quadratic form ||Z_g u||^2 = ||R u||^2, so
// the Profile keeps the direct-residual arithmetic (no normal-equation cancellation for the residuals).
// Trace terms (log det, sum c w, sum c^2 w^2) use the group sizes and power sums of delta.
// One compression of [UX0, UG_s, I_s, y] per SNP serves the genotype and interaction models (nested column
// subsets), and the grid of all models is one matrix product (grid coefficients x pair products).
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <array>
#include <limits>
using Eigen::ArrayXd; using Eigen::ArrayXXd; using Eigen::MatrixXd; using Eigen::VectorXd; using Eigen::VectorXi;
static const double NEG_INF = -std::numeric_limits<double>::infinity();
static const double L10 = 2.302585092994046;
static std::vector<double> make_grid() {
  const int nx = 80; const double xl = 5;
  std::vector<double> h(1, 0.0);
  for (int i = 0; i < nx; ++i) h.push_back(1 / (1 + std::pow(10.0, -(-xl + 2 * xl * i / (nx - 1.0)))));
  h.push_back(1); return h;
}
static const std::vector<double> HS = make_grid();
static const int NG = HS.size();
static const int MAX_PEAKS = 3;
static const int KN = 2;              // Taylor order of the null group
static double GROUP_TOL = 1e-12;      // set to -1 to disable grouping (all singletons: p1-like)
static int STATS[4];
#include <chrono>
static double TT[8];
static double now(){return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();}
// [[Rcpp::export]]
Rcpp::NumericVector timers(bool reset){Rcpp::NumericVector v(TT,TT+8); if(reset) for(auto&x:TT)x=0; return v;}                  // compressed evals, exact fallbacks, grid fallbacks

// ---------------- eigenvalue groups
struct Groups {
  int n = 0, ng = 0;
  std::vector<int> perm, start, size, K;      // contrasts in group order
  std::vector<double> cen, dmax;
  std::vector<std::array<double, KN + 1>> pw; // power sums of delta^k
};
static Groups make_groups(const ArrayXd& lam) {
  Groups G; const int n = lam.size(); G.n = n;
  std::vector<int> o(n); for (int i = 0; i < n; ++i) o[i] = i;
  std::stable_sort(o.begin(), o.end(), [&](int a, int b) { return lam[a] < lam[b]; });
  const double mx = lam.maxCoeff();
  int i = 0;
  while (i < n) {
    int j = i + 1; int K = 0;
    double c = lam[o[i]];
    if (GROUP_TOL < 0) { /* singleton */ }
    else if (c <= 1e-12 * mx) { K = KN; while (j < n && lam[o[j]] <= 1e-12 * mx) ++j; }
    else while (j < n && lam[o[j]] - c <= GROUP_TOL * c) ++j;
    G.start.push_back(i); G.size.push_back(j - i); G.cen.push_back(c);
    if (j - i == 1) K = 0;
    G.K.push_back(K);
    std::array<double, KN + 1> p{}; double dm = 0;
    for (int t = i; t < j; ++t) { double dl = lam[o[t]] - c, q = 1; dm = std::max(dm, dl); for (int k = 0; k <= KN; ++k) { p[k] += q; q *= dl; } }
    G.pw.push_back(p); G.dmax.push_back(dm);
    i = j;
  }
  G.perm = o; G.ng = G.start.size();
  return G;
}

// ---------------- compressed rows of Z (n x c)
struct Layout { std::vector<int> grp, ord; int m = 0; };   // row types, same for every Z with c columns
static Layout make_layout(const Groups& G, int c) {
  Layout L;
  for (int g = 0; g < G.ng; ++g)
    for (int k = 0; k <= G.K[g]; ++k) {
      if (k > 0 && G.dmax[g] == 0) break;
      int r = std::min(G.size[g], c);
      for (int t = 0; t < r; ++t) { L.grp.push_back(g); L.ord.push_back(k); }
    }
  L.m = L.grp.size(); return L;
}
static void compress(const Groups& G, const ArrayXd& lam, const MatrixXd& Z, MatrixXd& R) {
  const int c = Z.cols(); int row = 0;
  R.resize(0, 0); std::vector<MatrixXd> parts; int m = 0;
  MatrixXd blk;
  // first pass to size R
  for (int g = 0; g < G.ng; ++g) for (int k = 0; k <= G.K[g]; ++k) { if (k > 0 && G.dmax[g] == 0) break; m += std::min(G.size[g], c); }
  R.setZero(m, c);
  for (int g = 0; g < G.ng; ++g) {
    const int s0 = G.start[g], ns = G.size[g];
    for (int k = 0; k <= G.K[g]; ++k) {
      if (k > 0 && G.dmax[g] == 0) break;
      if (ns == 1) { R.row(row++) = Z.row(G.perm[s0]); continue; }
      blk.resize(ns, c);
      for (int t = 0; t < ns; ++t) {
        int i = G.perm[s0 + t];
        double sc = k == 0 ? 1.0 : std::pow(lam[i] - G.cen[g], 0.5 * k);
        blk.row(t) = sc * Z.row(i);
      }
      const int r = std::min(ns, c);
      Eigen::HouseholderQR<MatrixXd> qr(blk);
      R.block(row, 0, r, c) = qr.matrixQR().topRows(r).triangularView<Eigen::Upper>();
      row += r;
    }
  }
}

// ---------------- weight coefficients at h
// d = a + h c (a = 1 - h). Taylor coefficients in delta of w = 1/d, w^2, w^3: (-h)^k/d^(k+1), (k+1)(-h)^k/d^(k+2),
// (k+1)(k+2)/2 (-h)^k/d^(k+3); c(lambda) = lambda - 1 = (c_g - 1) + delta shifts the order.
struct Coef { ArrayXd w, cw2, c2w3; double logdet = 0, tcw = 0, tc2w2 = 0; bool ok = true; };
static void coefs(const Groups& G, const Layout& L, double h, bool deriv, Coef& C) {
  const int ng = G.ng; const double a = 1 - h;
  static thread_local std::vector<std::array<double, KN + 1>> W, W2, W3;   // workspace (no allocation per call)
  W.resize(ng); W2.resize(ng); W3.resize(ng);
  C.logdet = C.tcw = C.tc2w2 = 0; C.ok = true;
  for (int g = 0; g < ng; ++g) {
    double d = a + h * G.cen[g];
    if (!(d > 0)) { C.ok = false; return; }
    int K = G.K[g];
    if (K > 0 && std::pow(h * G.dmax[g] / d, K + 1) > 1e-13) { C.ok = false; return; }
    double q = 1 / d, mh = -h;
    for (int k = 0; k <= K; ++k) {
      W[g][k] = q; W2[g][k] = (k + 1) * q / d; W3[g][k] = 0.5 * (k + 1) * (k + 2) * q / (d * d);
      q *= mh / d;
    }
    const auto& p = G.pw[g];
    C.logdet += p[0] * std::log(d);
    double t = h / d, tk = t;
    for (int k = 1; k <= K; ++k) { C.logdet += ((k & 1) ? 1 : -1) * tk / k * p[k]; tk *= t; }
    if (deriv) {
      double l1 = G.cen[g] - 1;
      for (int k = 0; k <= K; ++k) {
        double cw = l1 * W[g][k] + (k ? W[g][k - 1] : 0);
        double c2w2 = l1 * l1 * W2[g][k] + (k ? 2 * l1 * W2[g][k - 1] : 0) + (k > 1 ? W2[g][k - 2] : 0);
        C.tcw += cw * p[k]; C.tc2w2 += c2w2 * p[k];
      }
    }
  }
  C.w.resize(L.m); if (deriv) { C.cw2.resize(L.m); C.c2w3.resize(L.m); }
  for (int r = 0; r < L.m; ++r) {
    int g = L.grp[r], k = L.ord[r];
    C.w[r] = W[g][k];
    if (deriv) {
      double l1 = G.cen[g] - 1;
      C.cw2[r] = l1 * W2[g][k] + (k ? W2[g][k - 1] : 0);
      C.c2w3[r] = l1 * l1 * W3[g][k] + (k ? 2 * l1 * W3[g][k - 1] : 0) + (k > 1 ? W3[g][k - 2] : 0);
    }
  }
}

// ---------------- Profile (baseline formulas), exact on full rows
struct Fit { double value = NAN, s2k = NAN, s2 = NAN; VectorXd b; bool ok = false; double se0 = NAN; };
struct Profile {
  ArrayXd d, w, c, r, v; MatrixXd Xw, A; VectorXd b, rhs, Xv; Eigen::LDLT<MatrixXd> L; double s = NAN;
  double operator()(const MatrixXd& X, const VectorXd& y, const ArrayXd& lam, double h, double* f1 = nullptr, double* f2 = nullptr) {
    const int n = y.size();
    d = h * lam + (1 - h);
    if ((d <= 0).any()) { if (f1) *f1 = *f2 = NAN; return NEG_INF; }
    w = d.inverse(); r = y.array();
    if (X.cols()) {
      Xw.noalias() = (X.array().colwise() * w).matrix();
      A.noalias() = Xw.transpose() * X; L.compute(A);
      rhs.noalias() = Xw.transpose() * y; b = L.solve(rhs); r -= (X * b).array();
    } else b.resize(0);
    ArrayXd Py = w * r; double g = (r * Py).sum(); s = g / n;
    if (!(s > 0)) { if (f1) *f1 = *f2 = NAN; return std::numeric_limits<double>::infinity(); }
    if (f1) {
      c = lam - 1; v = c * Py;
      double Q1 = (c * Py.square()).sum(), Q2 = (w * v.square()).sum();
      if (X.cols()) { Xv.noalias() = Xw.transpose() * v.matrix(); Q2 -= Xv.dot(L.solve(Xv)); }
      *f1 = 0.5 * n * Q1 / g - 0.5 * (c * w).sum();
      *f2 = -0.5 * n * (2 * Q2 / g - Q1 * Q1 / (g * g)) + 0.5 * (c.square() * w.square()).sum();
    }
    return -0.5 * n - 0.5 * (n * std::log(s) + d.log().sum());
  }
};

// Same on compressed rows (Xc, yc) with row coefficients; n = number of contrasts.
struct CProfile {
  const Groups* G; const Layout* Lay; Coef C; ArrayXd r, r2; MatrixXd Xw, A; VectorXd b, rhs, Xv; Eigen::LDLT<MatrixXd> L; double s = NAN;
  bool operator()(const MatrixXd& X, const VectorXd& y, double h, double& val, double* f1 = nullptr, double* f2 = nullptr) {
    const int n = G->n;
    coefs(*G, *Lay, h, f1 != nullptr, C);
    if (!C.ok) return false;
    const ArrayXd& w = C.w;
    r = y.array();
    if (X.cols()) {
      Xw.noalias() = (X.array().colwise() * w).matrix();
      A.noalias() = Xw.transpose() * X; L.compute(A);
      rhs.noalias() = Xw.transpose() * y; b = L.solve(rhs); r -= (X * b).array();
    } else b.resize(0);
    r2 = r.square();
    double g = (w * r2).sum(); s = g / n;
    if (!(s > 0)) { if (f1) *f1 = *f2 = NAN; val = std::numeric_limits<double>::infinity(); return true; }
    if (f1) {
      double Q1 = (C.cw2 * r2).sum(), Q2 = (C.c2w3 * r2).sum();
      if (X.cols()) { Xv.noalias() = X.transpose() * (C.cw2 * r).matrix(); Q2 -= Xv.dot(L.solve(Xv)); }
      *f1 = 0.5 * n * Q1 / g - 0.5 * C.tcw;
      *f2 = -0.5 * n * (2 * Q2 / g - Q1 * Q1 / (g * g)) + 0.5 * C.tc2w2;
    }
    val = -0.5 * n - 0.5 * (n * std::log(s) + C.logdet);
    return true;
  }
};

// A model: compressed rows and (lazily) the full rows for the exact fallback.
struct Model {
  const Groups* G; const Layout* Lay; const ArrayXd* lam;
  MatrixXd Xc; VectorXd yc;                 // compressed
  std::function<void(MatrixXd&, VectorXd&)> full; MatrixXd Xf; VectorXd yf; bool have_full = false;
  CProfile CP; Profile P; double s = NAN; VectorXd b; double ainv00 = NAN; int p0 = 0;
  double eval(double h, double* f1 = nullptr, double* f2 = nullptr, bool keep = false) {
    double v;
    CP.G = G; CP.Lay = Lay;
    if (CP(Xc, yc, h, v, f1, f2)) {
      ++STATS[0];
      if (keep) { s = CP.s; b = CP.b; if (Xc.cols() > p0) { VectorXd e = VectorXd::Unit(Xc.cols(), p0); ainv00 = CP.L.solve(e)[p0]; } }
      return v;
    }
    ++STATS[1];
    if (!have_full) { full(Xf, yf); have_full = true; }
    v = P(Xf, yf, *lam, h, f1, f2);
    if (keep) { s = P.s; b = P.b; if (Xf.cols() > p0) { VectorXd e = VectorXd::Unit(Xf.cols(), p0); ainv00 = P.L.solve(e)[p0]; } }
    return v;
  }
};

static std::pair<double, double> newton(Model& M, const double* v, int best) {
  const double XCAP = 12;
  auto X_of = [&](double h) { return h <= 0 ? -XCAP : h >= 1 ? XCAP : std::max(-XCAP, std::min(XCAP, std::log10(h / (1 - h)))); };
  auto H = [](double x) { return 1 / (1 + std::pow(10.0, -x)); };
  double hb = HS[best], vb = v[best], f1, f2;
  double a = HS[std::max(0, best - 1)], c = HS[std::min(NG - 1, best + 1)];
  if (!std::isfinite(vb) || !(c > a)) return {hb, vb};
  double lo = X_of(a), hi = X_of(c), x = X_of(hb);
  for (int it = 0; it < 60; ++it) {
    double h = H(x), fv = M.eval(h, &f1, &f2);
    if (fv > vb) { vb = fv; hb = h; }
    double dh = L10 * h * (1 - h), d2h = L10 * (1 - 2 * h) * dh;
    double g1 = f1 * dh, g2 = f2 * dh * dh + f1 * d2h;
    if (!std::isfinite(g1)) break;
    if (g1 > 0) lo = std::max(lo, x); else hi = std::min(hi, x);
    double xn = g2 < 0 ? x - g1 / g2 : NAN;
    if (!std::isfinite(xn) || xn <= lo || xn >= hi) xn = 0.5 * (lo + hi);
    if (std::fabs(xn - x) < 1e-9 || hi - lo < 1e-10) {
      double hh = H(xn), fv2 = M.eval(hh);
      if (fv2 > vb) { vb = fv2; hb = hh; }
      break;
    }
    x = xn;
  }
  return {hb, vb};
}

static Fit refine(Model& M, const double* v, int best) {
  std::vector<int> peaks;
  for (int i = 0; i < NG; ++i)
    if (std::isfinite(v[i]) && (i == 0 || v[i] >= v[i - 1]) && (i == NG - 1 || v[i] >= v[i + 1])) peaks.push_back(i);
  std::sort(peaks.begin(), peaks.end(), [&](int p, int q) { return v[p] > v[q]; });
  if (peaks.empty()) peaks.push_back(best);
  double hb = HS[peaks[0]], vb = NEG_INF;
  for (size_t j = 0; j < peaks.size() && j < (size_t)MAX_PEAKS; ++j) {
    auto m = newton(M, v, peaks[j]);
    if (m.second > vb) { vb = m.second; hb = m.first; }
  }
  Fit f; f.value = M.eval(hb, nullptr, nullptr, true); f.s2k = M.s * hb; f.s2 = M.s * (1 - hb); f.b = M.b;
  f.ok = std::isfinite(f.value);
  f.se0 = std::sqrt(M.s * M.ainv00);
  return f;
}

static int argmax(const double* v) { int b = 0; for (int i = 1; i < NG; ++i) if (v[i] > v[b]) b = i; return b; }

// Grid: GW (NG x m) row coefficients (rows of invalid h flagged in gok); pair products of the compressed rows of
// every SNP: Prod (m x S*npair). Then per (h, SNP, model) the normal equations as in the baseline grid.
struct Grid { MatrixXd GW; ArrayXd logdet; std::vector<bool> ok; };
static Grid make_grid_coef(const Groups& G, const Layout& L) {
  Grid g; g.GW.resize(NG, L.m); g.logdet.resize(NG); g.ok.resize(NG);
  Coef C;
  for (int i = 0; i < NG; ++i) {
    coefs(G, L, HS[i], false, C); g.ok[i] = C.ok;
    if (C.ok) { g.GW.row(i) = C.w.matrix().transpose(); g.logdet[i] = C.logdet; } else g.GW.row(i).setZero();
  }
  return g;
}
// Rc: list of compressed matrices (m x c), one per SNP. Returns sums SM (NG x S*npair), pair order upper-tri col-major.
static MatrixXd grid_sums(const Grid& gr, const std::vector<MatrixXd>& Rc) {
  const int S = Rc.size(); if (!S) return MatrixXd();
  const int m = Rc[0].rows(), c = Rc[0].cols(), np = c * (c + 1) / 2;
  MatrixXd Prod(m, (size_t)S * np);
  for (int s = 0; s < S; ++s) {
    int t = s * np;
    for (int j = 0; j < c; ++j) for (int l = 0; l <= j; ++l) Prod.col(t++) = Rc[s].col(j).cwiseProduct(Rc[s].col(l));
  }
  return gr.GW * Prod;
}
// value at grid point i for SNP s, model columns cols (indices < c-1), y = column c-1
static double grid_value(const MatrixXd& SM, int i, int s, int c, const std::vector<int>& cols, double logdet, int n) {
  const int np = c * (c + 1) / 2, p = cols.size();
  auto at = [&](int j, int l) { if (j < l) std::swap(j, l); return SM(i, (size_t)s * np + j * (j + 1) / 2 + l); };
  MatrixXd A(p, p); VectorXd rhs(p);
  for (int a = 0; a < p; ++a) { rhs[a] = at(cols[a], c - 1); for (int b = 0; b <= a; ++b) A(a, b) = A(b, a) = at(cols[a], cols[b]); }
  double yy = at(c - 1, c - 1);
  double rss = yy;
  if (p) { VectorXd b = A.ldlt().solve(rhs); rss -= rhs.dot(b); }
  double sv = rss / n;
  return sv > 0 ? -0.5 * n - 0.5 * (n * std::log(sv) + logdet) : std::numeric_limits<double>::infinity();
}

// Values of the nested models cols[0..len) for every len in lens, at grid point i, from one Cholesky of the
// augmented normal equations [A rhs; rhs' yy]: rss(len) = yy - sum_{j < len} z_j^2, z = L^-1 rhs.
static void grid_prefix(const MatrixXd& SM, int i, int s, int c, const std::vector<int>& cols, double logdet, int n,
                        const int* lens, int nl, double* out) {
  const int np = c * (c + 1) / 2, p = cols.size();
  const double* row = nullptr; (void)row;
  auto at = [&](int j, int l) { if (j < l) std::swap(j, l); return SM(i, (size_t)s * np + j * (j + 1) / 2 + l); };
  double Lm[16][16], z[16]; bool ok = p <= 16;
  double yy = at(c - 1, c - 1);
  for (int j = 0; ok && j < p; ++j) {
    for (int k = 0; k <= j; ++k) {
      double a = at(cols[j], cols[k]);
      for (int l = 0; l < k; ++l) a -= Lm[j][l] * Lm[k][l];
      if (k < j) Lm[j][k] = a / Lm[k][k];
      else if (a > 0) Lm[j][j] = std::sqrt(a); else ok = false;
    }
    if (!ok) break;
    double a = at(cols[j], c - 1);
    for (int l = 0; l < j; ++l) a -= Lm[j][l] * z[l];
    z[j] = a / Lm[j][j];
  }
  for (int t = 0; t < nl; ++t) {
    double rss = yy;
    if (ok) for (int j = 0; j < lens[t]; ++j) rss -= z[j] * z[j];
    else { std::vector<int> sub(cols.begin(), cols.begin() + lens[t]); out[t] = grid_value(SM, i, s, c, sub, logdet, n); continue; }
    double sv = rss / n;
    out[t] = sv > 0 ? -0.5 * n - 0.5 * (n * std::log(sv) + logdet) : std::numeric_limits<double>::infinity();
  }
}

// exact grid value (baseline formula) for the rare invalid points with all d > 0
static double grid_exact(const MatrixXd& X, const VectorXd& y, const ArrayXd& lam, double h) {
  ArrayXd d = h * lam + (1 - h); if ((d <= 0).any()) return NEG_INF;
  ArrayXd w = d.inverse(); const int n = y.size();
  double rss = (y.array().square() * w).sum();
  if (X.cols()) { MatrixXd Xw = (X.array().colwise() * w).matrix(); VectorXd rhs = Xw.transpose() * y; rss -= rhs.dot((Xw.transpose() * X).ldlt().solve(rhs)); }
  double sv = rss / n;
  return sv > 0 ? -0.5 * n - 0.5 * (n * std::log(sv) + d.log().sum()) : std::numeric_limits<double>::infinity();
}

static bool standardize(Eigen::Ref<VectorXd> x, double& sd) {
  const int n = x.size(); x.array() -= x.mean(); sd = std::sqrt(x.squaredNorm() / (n - 1));
  if (!(sd > 0)) return false;
  x /= sd; return true;
}

static Rcpp::List fit_gene(const MatrixXd& U, const VectorXd& lam_, const VectorXd& Uy, const MatrixXd& UX0,
                           const MatrixXd& G, const MatrixXd& C, bool boot) {
  const int S = G.cols(), q = C.cols(), p0 = UX0.cols(), n = Uy.size(); double t0 = now();
  const ArrayXd lam = lam_.array();
  const Groups GR = make_groups(lam);
  const bool anyzero = (lam <= 0).any();

  // null model on the data
  const int c0 = p0 + 1;
  Layout L0 = make_layout(GR, c0); Grid gr0 = make_grid_coef(GR, L0);
  MatrixXd Z0(n, c0); Z0 << UX0, Uy;
  Fit f0;
  {
    std::vector<MatrixXd> R(1); compress(GR, lam, Z0, R[0]);
    MatrixXd SM = grid_sums(gr0, R);
    std::vector<int> cols(p0); for (int j = 0; j < p0; ++j) cols[j] = j;
    VectorXd v(NG);
    for (int i = 0; i < NG; ++i) v[i] = gr0.ok[i] ? grid_value(SM, i, 0, c0, cols, gr0.logdet[i], n)
                                                : (anyzero ? NEG_INF : (++STATS[2], grid_exact(UX0, Uy, lam, HS[i])));
    Model M{&GR, &L0, &lam}; M.Xc = R[0].leftCols(p0); M.yc = R[0].col(p0); M.p0 = p0;
    M.full = [&](MatrixXd& X, VectorXd& y) { X = UX0; y = Uy; };
    f0 = refine(M, v.data(), argmax(v.data()));
  }
  if (!f0.ok) Rcpp::stop("no-genotype model did not converge");

  MatrixXd Gs = G; VectorXd gsd(S);
  for (int s = 0; s < S; ++s) standardize(Gs.col(s), gsd[s]);
  std::vector<MatrixXd> blocks(1, U.transpose() * Gs);
  std::vector<std::vector<bool>> iok(q, std::vector<bool>(S)); VectorXi df = VectorXi::Zero(S);
  for (int j = 0; j < q; ++j) {
    MatrixXd I = G.array().colwise() * C.col(j).array();
    for (int s = 0; s < S; ++s) { double sd; iok[j][s] = standardize(I.col(s), sd); if (!iok[j][s]) I.col(s).setZero(); else ++df[s]; }
    blocks.push_back(U.transpose() * I);
  }
  TT[0]+=now()-t0; t0=now();
  const MatrixXd& UG = blocks[0];
  // full column set per SNP: [UX0, UG_s, I_1s..I_qs, y]
  const int c = p0 + 1 + q + 1;
  auto Zfull = [&](int s, const VectorXd& y) { MatrixXd Z(n, c); Z.leftCols(p0) = UX0; Z.col(p0) = UG.col(s);
    for (int j = 0; j < q; ++j) Z.col(p0 + 1 + j) = blocks[1 + j].col(s); Z.col(c - 1) = y; return Z; };
  auto cols1 = [&]() { std::vector<int> k(p0 + 1); for (int j = 0; j <= p0; ++j) k[j] = j; return k; };
  auto cols2 = [&](int s) { auto k = cols1(); for (int j = 0; j < q; ++j) if (iok[j][s]) k.push_back(p0 + 1 + j); return k; };
  std::vector<int> colsX0(p0); for (int j = 0; j < p0; ++j) colsX0[j] = j;
  Layout L = make_layout(GR, c); Grid gr = make_grid_coef(GR, L);

  // grid values of the nested models (prefixes of cols2(s)) for SNP s: V (NG x 3) for lens {p0, p0+1, p0+1+df}
  auto grid3 = [&](const MatrixXd& SM, int s, const VectorXd& y, MatrixXd& V) {
    auto cols = cols2(s); int lens[3] = {p0, p0 + 1, (int)cols.size()};
    V.resize(NG, 3);
    for (int i = 0; i < NG; ++i) {
      double o[3];
      if (gr.ok[i]) grid_prefix(SM, i, s, c, cols, gr.logdet[i], n, lens, 3, o);
      else for (int t = 0; t < 3; ++t) {
        if (anyzero) { o[t] = NEG_INF; continue; }
        ++STATS[2]; MatrixXd Z = Zfull(s, y), X(n, lens[t]); for (int j = 0; j < lens[t]; ++j) X.col(j) = Z.col(cols[j]);
        o[t] = grid_exact(X, y, lam, HS[i]);
      }
      for (int t = 0; t < 3; ++t) V(i, t) = o[t];
    }
  };
  auto run_model = [&](const MatrixXd& R, const std::vector<int>& cols, const double* vg, int s, const VectorXd& y) {
    VectorXd v = Eigen::Map<const VectorXd>(vg, NG);
    const int p = cols.size();
    MatrixXd Xc(R.rows(), p); for (int j = 0; j < p; ++j) Xc.col(j) = R.col(cols[j]);
    auto fullf = [&, s](MatrixXd& X, VectorXd& yy) { MatrixXd Z = Zfull(s, y); X.resize(n, p); for (int j = 0; j < p; ++j) X.col(j) = Z.col(cols[j]); yy = y; };
    Model M{&GR, &L, &lam}; M.Xc = Xc; M.yc = R.col(c - 1); M.p0 = p0; M.full = fullf;
    return refine(M, v.data(), argmax(v.data()));
  };

  auto ok_nested = [](Fit& f, const Fit& g0) { if (!f.ok || f.value < g0.value - 1e-6) f.ok = false; };
  Rcpp::NumericMatrix out(S, 7);
  std::fill(out.begin(), out.end(), NA_REAL);

  std::vector<MatrixXd> R(S);
  t0=now();
  for (int s = 0; s < S; ++s) compress(GR, lam, Zfull(s, Uy), R[s]);
  TT[1]+=now()-t0; t0=now();
  MatrixXd SM = grid_sums(gr, R);
  TT[2]+=now()-t0; t0=now();
  std::vector<Fit> f1(S);
  for (int s = 0; s < S; ++s) {
    out(s, 2) = df[s];
    double tg=now(); MatrixXd V; grid3(SM, s, Uy, V); TT[5]+=now()-tg;
    f1[s] = run_model(R[s], cols1(), V.col(1).data(), s, Uy); ok_nested(f1[s], f0);
    if (!f1[s].ok) continue;
    out(s, 0) = f1[s].value;
    out(s, 3) = f1[s].b[p0] / gsd[s]; out(s, 4) = f1[s].se0 / gsd[s];
    if (df[s] > 0) {
      Fit f2 = run_model(R[s], cols2(s), V.col(2).data(), s, Uy); ok_nested(f2, f1[s]);
      if (f2.ok) out(s, 1) = f2.value;
    }
  }
  TT[3]+=now()-t0; t0=now();
  if (!boot) return Rcpp::List::create(Rcpp::Named("l0") = f0.value, Rcpp::Named("res") = out);

  Rcpp::RNGScope rng;
  MatrixXd Yb = MatrixXd::Zero(n, S);
  for (int s = 0; s < S; ++s) {
    if (!f1[s].ok) continue;
    for (int i = 0; i < n; ++i) Yb(i, s) = std::sqrt(f1[s].s2k * lam[i] + f1[s].s2) * R::norm_rand();
    Yb.col(s) += UX0 * f1[s].b.head(p0) + UG.col(s) * f1[s].b[p0];
  }
  TT[4]+=now()-t0; t0=now();
  for (int s = 0; s < S; ++s) compress(GR, lam, Zfull(s, Yb.col(s)), R[s]);
  TT[1]+=now()-t0; t0=now();
  SM = grid_sums(gr, R);
  TT[2]+=now()-t0; t0=now();
  for (int s = 0; s < S; ++s) {
    if (!f1[s].ok) continue;
    VectorXd y = Yb.col(s);
    double tg=now(); MatrixXd V; grid3(SM, s, y, V); TT[5]+=now()-tg;
    Fit g0 = run_model(R[s], colsX0, V.col(0).data(), s, y);
    if (!g0.ok) continue;
    Fit g1 = run_model(R[s], cols1(), V.col(1).data(), s, y); ok_nested(g1, g0);
    if (!g1.ok) continue;
    out(s, 5) = g1.value;
    if (df[s] > 0) {
      Fit g2 = run_model(R[s], cols2(s), V.col(2).data(), s, y); ok_nested(g2, g1);
      if (g2.ok) out(s, 6) = g2.value;
    }
  }
  TT[6]+=now()-t0;
  return Rcpp::List::create(Rcpp::Named("l0") = f0.value, Rcpp::Named("res") = out);
}

// [[Rcpp::export]]
Rcpp::List fit_gene_proto(Eigen::MatrixXd U, Eigen::VectorXd lambda, Eigen::VectorXd Uy, Eigen::MatrixXd UX0,
                          Eigen::MatrixXd G, Eigen::MatrixXd C, bool boot) { return fit_gene(U, lambda, Uy, UX0, G, C, boot); }
// [[Rcpp::export]]
Rcpp::IntegerVector proto_stats(bool reset) { Rcpp::IntegerVector v(STATS, STATS + 3); if (reset) STATS[0] = STATS[1] = STATS[2] = 0; return v; }
// [[Rcpp::export]]
void proto_group_tol(double t) { GROUP_TOL = t; }
