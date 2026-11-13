// [[Rcpp::depends(RcppEigen)]]
// p9: combined numerical candidate (p7 + reviewer round-3 changes). Same model, grid, peaks and Newton steps as
// src/lmm.cpp; only how the weighted sums are formed changes.
//
// Model on rotated data: Uy = UX b + e, e ~ N(0, diag(s (h lam + (1 - h)))). Every quantity the fits need is a
// weighted sum over contrasts t, sum_t f(h, lam_t) u_t v_t, with f = w = 1 / d, c w^2 or c^2 w^3 (c = lam - 1,
// d = h lam + (1 - h)). Contrasts are grouped by eigenvalue ("rows" below) so these become sum_r f_r(h) M_r(u, v)
// with moments M_r computed once per gene / SNP:
//  * non-null eigenvalues tied within a relative 1e-12 share one weight (taken at the group mean);
//  * the null group (eigen-decomposition noise around 0) is kept exact: f(c + delta) is expanded in delta to
//    order KN, one row per power (moments sum delta^k u v); where the truncation is not negligible at some h,
//    its contrasts are summed one by one at that h.
// The grid evaluates all h at once (moments x grid weights, a GEMM) and solves the per-SNP normal equations by
// the Schur complement on the shared X0 block. Where the normal equations cancel, values are recomputed from
// direct residuals. SNPs are processed in blocks (rotation, moments, grid, fits, bootstrap), so memory does not
// grow with the moments of a whole gene. Bootstrap normals are drawn in SNP and contrast order, as in p0.
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
using Eigen::ArrayXd; using Eigen::MatrixXd; using Eigen::VectorXd; using Eigen::VectorXi;
static const double NEG_INF = -std::numeric_limits<double>::infinity();
static const double POS_INF = std::numeric_limits<double>::infinity();
static const double L10 = 2.302585092994046;   // log(10)

// h grid as p0: h = 0, x = log10(h / (1 - h)) uniform on [-5, 5] (80 points), h = 1.
static std::vector<double> make_grid() {
  const int nx = 80; const double xl = 5;
  std::vector<double> h(1, 0.0);
  for (int i = 0; i < nx; ++i) h.push_back(1 / (1 + std::pow(10.0, -(-xl + 2 * xl * i / (nx - 1.0)))));
  h.push_back(1); return h;
}
static const std::vector<double> HS = make_grid();
static const int NG = HS.size();
static const int MAX_PEAKS = 3;   // local maxima of the grid refined per fit (as p0)

// Thresholds (why each value):
// Non-null eigenvalues within this relative spread share one weight: |w(lam) - w(lam_g)| / w <= h |dlam| / d
// <= |dlam| / lam <= 1e-12 at every h, so the log-likelihood moves by <= n 1e-12 (below the 1e-9 target).
static const double GROUP_RTOL = 1e-12;
// Eigenvalues <= this * max are numerical zeros of the kernel (R's eigen() leaves values up to ~1e-14 max):
// they are the null group, kept exact by the expansion below rather than tied.
static const double NULL_RTOL = 1e-12;
// Order of the expansion of f(c + delta) in delta for the null group, and the bound on its truncation:
// with rho = h max(delta) / d the relative error of every null weight is <= rho^(KN+1) <= 1e-13 (order 3 keeps
// the expansion valid up to x ~ 11 on the default kernel and x ~ 5.4 on a kernel scaled by 1e6).
static const int KN = 3;
static const double TAYLOR_TOL = 1e-13;
// Grid: rss = y'Wy - rhs'b loses all digits when rss / y'Wy approaches eps; below 1e-10 the grid value is
// recomputed from direct residuals (reviewer: near-perfect fits).
static const double CANCEL_GRID = 1e-10;
// Profile (Newton values, final values): g = coef . moments has relative error ~ eps sum|terms| / g; below 1e-3
// (error >~ 2e-13) it switches to direct residuals, as p0 always does.
static const double CANCEL_PROF = 2;
// Schur pivot (relative to the column's own weighted norm) below which the per-SNP system is solved by p0's
// LDLT on the full matrix (collinear SNP / interaction / covariates: p0's pseudo-solve).
static const double SCHUR_PIVOT = 1e-10;
// SNPs per block: bounds the moments (rows x pairs x block) and the rotated block, with GEMMs still wide.
static const int SNP_BLOCK = 128;

static inline int pidx(int a, int b, int N) { if (a > b) std::swap(a, b); return a * N - a * (a - 1) / 2 + (b - a); }

// ---------------------------------------------------------------- eigenvalue groups
struct Weights {   // per-row coefficients of w, c w^2, c^2 w^3 at one h, and the trace terms
  ArrayXd a, b, e; double T1 = 0, T2 = 0, LD = 0; bool null_ok = true;
};

struct Spectrum {
  struct Row { int t0, t1, k; double lam; };      // contrasts [t0, t1) (permuted order); k = power of delta
  int n, n0 = 0, K = 0; double lmin, lmax, c = 0, dmax = 0;
  std::vector<int> perm; std::vector<Row> rows;
  std::vector<std::vector<double>> dpow;         // dpow[k][t] = delta_t^k, t < n0 (null contrasts come first)
  std::vector<double> lamnull; double P[KN + 1] = {0};   // exact null eigenvalues; power sums of delta

  explicit Spectrum(const ArrayXd& lam) : n(lam.size()), perm(lam.size()) {
    for (int i = 0; i < n; ++i) perm[i] = i;
    std::stable_sort(perm.begin(), perm.end(), [&](int a, int b) { return lam[a] < lam[b]; });
    lmin = lam.minCoeff(); lmax = lam.maxCoeff();
    while (n0 < n && lam[perm[n0]] <= NULL_RTOL * std::max(lmax, 0.0)) ++n0;
    if (n0) {
      c = lam[perm[0]];
      for (int t = 0; t < n0; ++t) { lamnull.push_back(lam[perm[t]]); dmax = std::max(dmax, lam[perm[t]] - c); }
      K = dmax > 0 ? KN : 0;
      dpow.assign(K + 1, std::vector<double>(n0, 1.0));
      for (int k = 1; k <= K; ++k) for (int t = 0; t < n0; ++t) dpow[k][t] = dpow[k - 1][t] * (lamnull[t] - c);
      for (int k = 0; k <= K; ++k) { for (int t = 0; t < n0; ++t) P[k] += dpow[k][t]; rows.push_back({0, n0, k, c}); }
    }
    for (int t = n0; t < n;) {
      const double l0 = lam[perm[t]]; double sum = 0; int u = t;
      while (u < n && lam[perm[u]] - l0 <= GROUP_RTOL * l0) sum += lam[perm[u++]];
      rows.push_back({t, u, 0, sum / (u - t)}); t = u;
    }
  }
  int R() const { return rows.size(); }
  int nnull_rows() const { return n0 ? K + 1 : 0; }
  // p0's any(d <= 0): d is linear in lam, so its extremes decide
  bool is_bad(double h) const { return !(h * lmin + (1 - h) > 0 && h * lmax + (1 - h) > 0); }
  MatrixXd permute(const MatrixXd& M) const { MatrixXd P(M.rows(), M.cols()); for (int i = 0; i < n; ++i) P.row(i) = M.row(perm[i]); return P; }
  // moments of the pair (u, v) on every row
  void mom(const double* u, const double* v, double* out) const {
    for (int r = 0; r < R(); ++r) {
      const Row& q = rows[r]; double acc = 0;
      if (q.k == 0) for (int t = q.t0; t < q.t1; ++t) acc += u[t] * v[t];
      else { const double* dk = dpow[q.k].data(); for (int t = q.t0; t < q.t1; ++t) acc += dk[t] * u[t] * v[t]; }
      out[r] = acc;
    }
  }
  // null contrast t's exact weights at h
  void null_w(int t, double h, double& w, double& wb, double& we) const {
    const double cl = lamnull[t] - 1; w = 1 / (h * lamnull[t] + (1 - h)); wb = cl * w * w; we = cl * wb * cl * w;
  }
  void weights(double h, bool derivs, Weights& W) const {
    const int R_ = R(), nr = nnull_rows();
    W.a.setZero(R_); W.b.setZero(R_); W.e.setZero(R_); W.T1 = W.T2 = W.LD = 0; W.null_ok = true;
    for (int r = nr; r < R_; ++r) {
      const double lg = rows[r].lam, d = h * lg + (1 - h), w = 1 / d, cg = lg - 1, cnt = rows[r].t1 - rows[r].t0;
      W.a[r] = w; W.LD += cnt * std::log(d);
      if (derivs) { W.b[r] = cg * w * w; W.e[r] = cg * cg * w * w * w; W.T1 += cnt * cg * w; W.T2 += cnt * cg * w * cg * w; }
    }
    if (!n0) return;
    const double d = h * c + (1 - h);
    W.null_ok = dmax == 0 || std::pow(h * dmax / d, KN + 1) <= TAYLOR_TOL;
    if (!W.null_ok) {   // sum the null contrasts one by one (callers add their pair sums directly)
      for (int t = 0; t < n0; ++t) {
        double w, wb, we; null_w(t, h, w, wb, we); const double cl = lamnull[t] - 1;
        W.LD += std::log(h * lamnull[t] + (1 - h)); W.T1 += cl * w; W.T2 += cl * w * cl * w;
      }
      return;
    }
    // f(c + delta) = sum_k f_k delta^k: w^m = d^-m sum_k C(m+k-1, k) q^k (q = -h delta / d); c(lam) = c0 + delta
    double W1[KN + 3] = {0}, W2[KN + 3] = {0}, W3[KN + 3] = {0}; const double q = -h / d, c0 = c - 1;
    double qk = 1;
    for (int k = 0; k <= K; ++k, qk *= q) { W1[k + 2] = qk / d; W2[k + 2] = (k + 1) * qk / (d * d); W3[k + 2] = 0.5 * (k + 1) * (k + 2) * qk / (d * d * d); }
    double hd = 1;
    for (int k = 0; k <= K; ++k) {   // index k + 2: W*[k + 1] and W*[k] are the k-1 and k-2 terms (0 below 0)
      const int i = k + 2;
      W.a[k] = W1[i]; W.LD += (k == 0 ? std::log(d) : ((k % 2) ? 1 : -1) * hd / k) * P[k];
      if (derivs) {
        W.b[k] = c0 * W2[i] + W2[i - 1]; W.e[k] = c0 * c0 * W3[i] + 2 * c0 * W3[i - 1] + W3[i - 2];
        W.T1 += (c0 * W1[i] + W1[i - 1]) * P[k]; W.T2 += (c0 * c0 * W2[i] + 2 * c0 * W2[i - 1] + W2[i - 2]) * P[k];
      }
      hd *= h / d;
    }
  }
  // sum_t f_j(h, lam_t) u_t v_t for the weight set j (0: w, 1: c w^2, 2: c^2 w^3), on contrasts directly
  double wsum(const Weights& W, int j, double h, const double* u, const double* v) const {
    const ArrayXd& f = j == 0 ? W.a : j == 1 ? W.b : W.e; double acc = 0;
    for (int r = 0; r < R(); ++r) {
      if (f[r] == 0) continue;
      const Row& q = rows[r]; double s = 0;
      if (q.k == 0) for (int t = q.t0; t < q.t1; ++t) s += u[t] * v[t];
      else { const double* dk = dpow[q.k].data(); for (int t = q.t0; t < q.t1; ++t) s += dk[t] * u[t] * v[t]; }
      acc += f[r] * s;
    }
    if (!W.null_ok) for (int t = 0; t < n0; ++t) { double w[3]; null_w(t, h, w[0], w[1], w[2]); acc += w[j] * u[t] * v[t]; }
    return acc;
  }
};

// Grid weights of a gene: row coefficients (NG x R) of w, with the null group's direct flag, and the exact
// log determinant per grid point.
struct GridW {
  MatrixXd A; std::vector<Weights> W; std::vector<char> bad; VectorXd logdet;
  GridW(const Spectrum& sp, const ArrayXd& lam) : A(NG, sp.R()), W(NG), bad(NG), logdet(NG) {
    for (int i = 0; i < NG; ++i) {
      bad[i] = sp.is_bad(HS[i]);
      if (bad[i]) { A.row(i).setZero(); continue; }
      sp.weights(HS[i], false, W[i]); A.row(i) = W[i].a.matrix().transpose();
      logdet[i] = (HS[i] * lam + (1 - HS[i])).log().sum();   // as p0
    }
  }
};

// ---------------------------------------------------------------- moments of one SNP block
// Columns [X0, B_1..B_k, Y]; X0 x X0 pair moments come from the gene (MX0), pairs involving only X0 and a shared
// y once per block (MY), every other pair per SNP (Mper, rows x nper*S).
struct Block {
  const Spectrum& sp; const MatrixXd& X0; const MatrixXd& MX0; const std::vector<MatrixXd>& B; const MatrixXd& Y;
  int p0, k, N, iy, S, nper = 0; bool shared_y;
  std::vector<int> kind, slot;   // per pair: 0 = MX0, 1 = MY, 2 = Mper
  MatrixXd MY, Mper;
  const double* col(int a, int s) const { return a < p0 ? X0.col(a).data() : a < iy ? B[a - p0].col(s).data() : Y.col(shared_y ? 0 : s).data(); }
  Block(const Spectrum& sp_, const MatrixXd& X0_, const MatrixXd& MX0_, const std::vector<MatrixXd>& B_, const MatrixXd& Y_, int S_)
      : sp(sp_), X0(X0_), MX0(MX0_), B(B_), Y(Y_), p0(X0_.cols()), k(B_.size()), N(p0 + k + 1), iy(p0 + k), S(S_), shared_y(Y_.cols() == 1) {
    const int np = N * (N + 1) / 2; kind.resize(np); slot.resize(np); int ny = 0;
    auto shared = [&](int a) { return a < p0 || (a == iy && shared_y); };
    for (int a = 0; a < N; ++a) for (int b = a; b < N; ++b) {
      const int t = pidx(a, b, N);
      if (b < p0) { kind[t] = 0; slot[t] = pidx(a, b, p0); }
      else if (shared(a) && shared(b)) { kind[t] = 1; slot[t] = ny++; }
      else { kind[t] = 2; slot[t] = nper++; }
    }
    MY.resize(sp.R(), ny); Mper.resize(sp.R(), (size_t)nper * S);
    for (int a = 0; a < N; ++a) for (int b = a; b < N; ++b) {
      const int t = pidx(a, b, N);
      if (kind[t] == 1) sp.mom(col(a, 0), col(b, 0), MY.col(slot[t]).data());
      else if (kind[t] == 2) for (int s = 0; s < S; ++s) sp.mom(col(a, s), col(b, s), Mper.col((size_t)s * nper + slot[t]).data());
    }
  }
  // moment column (rows) of pair (a, b) for SNP s
  const double* pm(int a, int b, int s) const {
    const int t = pidx(a, b, N);
    return kind[t] == 0 ? MX0.col(slot[t]).data() : kind[t] == 1 ? MY.col(slot[t]).data() : Mper.col((size_t)s * nper + slot[t]).data();
  }
};

static MatrixXd x0_moments(const Spectrum& sp, const MatrixXd& X0) {
  const int p0 = X0.cols(); MatrixXd M(sp.R(), p0 * (p0 + 1) / 2);
  for (int a = 0; a < p0; ++a) for (int b = a; b < p0; ++b) sp.mom(X0.col(a).data(), X0.col(b).data(), M.col(pidx(a, b, p0)).data());
  return M;
}

// A model of SNP s: columns cols (indices into the block's columns) and y, as moments (rows x pairs of
// [cols, y]) plus the raw columns for the direct-residual and direct-null paths.
struct Model { MatrixXd MM; std::vector<const double*> z; };   // z: columns then y
static Model make_model(const Block& bl, int s, const std::vector<int>& cols) {
  const int p = cols.size(), Nm = p + 1; Model m; m.MM.resize(bl.sp.R(), Nm * (Nm + 1) / 2);
  auto c = [&](int a) { return a < p ? cols[a] : bl.iy; };
  for (int a = 0; a < Nm; ++a) {
    m.z.push_back(bl.col(c(a), s));
    for (int b = a; b < Nm; ++b) m.MM.col(pidx(a, b, Nm)) = Eigen::Map<const VectorXd>(bl.pm(c(a), c(b), s), bl.sp.R());
  }
  return m;
}

// ---------------------------------------------------------------- profile likelihood (Newton, final values)
// As p0: f = -(n/2) log(g/n) - (1/2) sum log d - n/2, with g = r'Wr, f' and f'' from Q1 = sum c w^2 r^2,
// Q2 = sum c^2 w^3 r^2 - Xv'A^-1 Xv, Xv = X' diag(c w^2) r, the trace terms sum c w and sum c^2 w^2.
// Sums come from the pair moments (MM' x row weights); r'Dr = coef . (MM' D) with coef the quadratic form of
// r = y - Xb in the pairs; with direct residuals where that cancels.
struct Profile {
  const Spectrum& sp; const int n; Weights W; MatrixXd AB, MW, A; VectorXd b, rhs, coef, Xv, tmp; ArrayXd r;
  Eigen::LDLT<MatrixXd> L; double s = NAN;
  explicit Profile(const Spectrum& s_) : sp(s_), n(s_.n) {}
  double operator()(const Model& m, double h, double* f1 = nullptr, double* f2 = nullptr) {
    const MatrixXd& MM = m.MM; const int Nm = m.z.size(), p = Nm - 1, nc = f1 ? 3 : 1, P2 = MM.cols();
    if (sp.is_bad(h)) { if (f1) *f1 = *f2 = NAN; return NEG_INF; }
    sp.weights(h, f1 != nullptr, W);
    AB.resize(sp.R(), 3); MW.resize(P2, 3);   // fixed shapes: no reallocation between value-only and derivative calls
    AB.col(0) = W.a.matrix(); if (f1) { AB.col(1) = W.b.matrix(); AB.col(2) = W.e.matrix(); }
    MW.leftCols(nc).noalias() = MM.transpose() * AB.leftCols(nc);
    if (!W.null_ok)
      for (int t = 0; t < sp.n0; ++t) {
        double w[3]; sp.null_w(t, h, w[0], w[1], w[2]);
        for (int a = 0; a < Nm; ++a) for (int e = a; e < Nm; ++e) {
          const double zz = m.z[a][t] * m.z[e][t]; const int q = pidx(a, e, Nm);
          for (int j = 0; j < nc; ++j) MW(q, j) += w[j] * zz;
        }
      }
    coef.setZero(P2); coef[pidx(p, p, Nm)] = 1;
    if (p) {
      A.resize(p, p); rhs.resize(p);
      for (int a = 0; a < p; ++a) { rhs[a] = MW(pidx(a, p, Nm), 0); for (int e = 0; e <= a; ++e) A(a, e) = A(e, a) = MW(pidx(e, a, Nm), 0); }
      L.compute(A); b = L.solve(rhs);
      for (int a = 0; a < p; ++a) { coef[pidx(a, p, Nm)] = -2 * b[a]; for (int e = a; e < p; ++e) coef[pidx(a, e, Nm)] = (a == e ? 1 : 2) * b[a] * b[e]; }
    } else b.resize(0);
    double g = coef.dot(MW.col(0));
    // the terms of g can be large with opposite signs: near-perfect fits, or b huge along a null direction of a
    // rank-deficient X (LDLT pseudo-solve); measure the cancellation against the sum of |terms|
    const bool direct = !(g > CANCEL_PROF * coef.cwiseAbs().dot(MW.col(0).cwiseAbs()));
    if (direct) {
      r = Eigen::Map<const ArrayXd>(m.z[p], n);
      for (int a = 0; a < p; ++a) r -= b[a] * Eigen::Map<const ArrayXd>(m.z[a], n);
      g = sp.wsum(W, 0, h, r.data(), r.data());
    }
    s = g / n;
    if (!(s > 0)) { if (f1) *f1 = *f2 = NAN; return POS_INF; }
    if (f1) {
      double Q1, Q2; Xv.resize(p);
      if (direct) {
        Q1 = sp.wsum(W, 1, h, r.data(), r.data()); Q2 = sp.wsum(W, 2, h, r.data(), r.data());
        for (int a = 0; a < p; ++a) Xv[a] = sp.wsum(W, 1, h, m.z[a], r.data());
      } else {
        Q1 = coef.dot(MW.col(1)); Q2 = coef.dot(MW.col(2));
        for (int a = 0; a < p; ++a) { double x = MW(pidx(a, p, Nm), 1); for (int e = 0; e < p; ++e) x -= MW(pidx(a, e, Nm), 1) * b[e]; Xv[a] = x; }
      }
      if (p) { tmp = L.solve(Xv); Q2 -= Xv.dot(tmp); }
      *f1 = 0.5 * n * Q1 / g - 0.5 * W.T1;
      *f2 = -0.5 * n * (2 * Q2 / g - Q1 * Q1 / (g * g)) + 0.5 * W.T2;
    }
    return -0.5 * n - 0.5 * (n * std::log(s) + W.LD);
  }
};

struct Fit { double value = NAN, s2k = NAN, s2 = NAN; VectorXd b; bool ok = false; };

// Safeguarded Newton in x = log10(delta) around grid point `best`, as p0.
static std::pair<double, double> newton(Profile& P, const Model& m, const double* v, int best) {
  const double XCAP = 12;
  auto X_of = [&](double h) { return h <= 0 ? -XCAP : h >= 1 ? XCAP : std::max(-XCAP, std::min(XCAP, std::log10(h / (1 - h)))); };
  auto H = [](double x) { return 1 / (1 + std::pow(10.0, -x)); };
  double hb = HS[best], vb = v[best], f1, f2;
  double a = HS[std::max(0, best - 1)], c = HS[std::min(NG - 1, best + 1)];
  if (!std::isfinite(vb) || !(c > a)) return {hb, vb};
  double lo = X_of(a), hi = X_of(c), x = X_of(hb);
  for (int it = 0; it < 60; ++it) {
    double h = H(x), fv = P(m, h, &f1, &f2);
    if (fv > vb) { vb = fv; hb = h; }
    double dh = L10 * h * (1 - h), d2h = L10 * (1 - 2 * h) * dh;
    double g1 = f1 * dh, g2 = f2 * dh * dh + f1 * d2h;
    if (!std::isfinite(g1)) break;
    if (g1 > 0) lo = std::max(lo, x); else hi = std::min(hi, x);
    double xn = g2 < 0 ? x - g1 / g2 : NAN;
    if (!std::isfinite(xn) || xn <= lo || xn >= hi) xn = 0.5 * (lo + hi);
    if (std::fabs(xn - x) < 1e-9 || hi - lo < 1e-10) {
      double hh = H(xn), fv2 = P(m, hh);
      if (fv2 > vb) { vb = fv2; hb = hh; }
      break;
    }
    x = xn;
  }
  return {hb, vb};
}

static int argmax(const double* v) { int b = 0; for (int i = 1; i < NG; ++i) if (v[i] > v[b]) b = i; return b; }

// Refine every local maximum of the grid values v (up to MAX_PEAKS, highest first); P is left at the best fit.
static Fit refine(Profile& P, const Model& m, const double* v) {
  std::vector<int> peaks;
  for (int i = 0; i < NG; ++i)
    if (std::isfinite(v[i]) && (i == 0 || v[i] >= v[i - 1]) && (i == NG - 1 || v[i] >= v[i + 1])) peaks.push_back(i);
  std::sort(peaks.begin(), peaks.end(), [&](int p, int q) { return v[p] > v[q]; });
  if (peaks.empty()) peaks.push_back(argmax(v));
  double hb = HS[peaks[0]], vb = NEG_INF;
  for (size_t j = 0; j < peaks.size() && j < (size_t)MAX_PEAKS; ++j) {
    auto r = newton(P, m, v, peaks[j]);
    if (r.second > vb) { vb = r.second; hb = r.first; }
  }
  Fit f; f.value = P(m, hb); f.s2k = P.s * hb; f.s2 = P.s * (1 - hb); f.b = P.b;
  f.ok = std::isfinite(f.value);
  return f;
}

// ---------------------------------------------------------------- grid values
// V (S x NG) for each model of `models` (block-column indices added to X0; a zero column is left out for that
// SNP as in p0). R = moments x grid weights; the null group's direct sums are added at grid points where its
// expansion does not hold. Normal equations by the Schur complement on X0 (A00^-1 per grid point, shared),
// the k x k remainder by Cholesky; LDLT on the full system (p0) for a small pivot; direct residuals where
// rss / y'Wy < CANCEL_GRID.
struct X0Grid { std::vector<MatrixXd> A0i; };
static X0Grid x0_grid(const Spectrum& sp, const GridW& gw, const MatrixXd& X0, const MatrixXd& MX0) {
  const int p0 = X0.cols(); X0Grid xg; xg.A0i.resize(NG);
  const MatrixXd R = MX0.transpose() * gw.A.transpose();
  for (int i = 0; i < NG; ++i) {
    if (gw.bad[i]) continue;
    MatrixXd A00(p0, p0);
    for (int a = 0; a < p0; ++a) for (int b = 0; b < p0; ++b) A00(a, b) = gw.W[i].null_ok ? R(pidx(a, b, p0), i) : sp.wsum(gw.W[i], 0, HS[i], X0.col(a).data(), X0.col(b).data());
    xg.A0i[i] = A00.ldlt().solve(MatrixXd::Identity(p0, p0));
  }
  return xg;
}

static double grid_direct(const Block& bl, const GridW& gw, int i, int s, const std::vector<int>& cols) {
  const int p = cols.size(), n = bl.sp.n; const double h = HS[i]; const Weights& W = gw.W[i];
  std::vector<const double*> z; for (int c : cols) z.push_back(bl.col(c, s));
  const double* y = bl.col(bl.iy, s);
  MatrixXd A(p, p); VectorXd rhs(p);
  for (int a = 0; a < p; ++a) { rhs[a] = bl.sp.wsum(W, 0, h, z[a], y); for (int e = 0; e <= a; ++e) A(a, e) = A(e, a) = bl.sp.wsum(W, 0, h, z[a], z[e]); }
  for (int a = 0; a < p; ++a) if (A(a, a) == 0) { A(a, a) = 1; rhs[a] = 0; }
  VectorXd b = A.ldlt().solve(rhs); ArrayXd r = Eigen::Map<const ArrayXd>(y, n);
  for (int a = 0; a < p; ++a) r -= b[a] * Eigen::Map<const ArrayXd>(z[a], n);
  return bl.sp.wsum(W, 0, h, r.data(), r.data());
}

static std::vector<MatrixXd> grid_values(const Block& bl, const GridW& gw, const X0Grid& xg, const std::vector<std::vector<int>>& models) {
  const Spectrum& sp = bl.sp; const int p0 = bl.p0, k = bl.k, iy = bl.iy, S = bl.S, n = sp.n, ny = bl.MY.cols(), nper = bl.nper;
  const MatrixXd At = gw.A.transpose();
  MatrixXd RY = bl.MY.transpose() * At, RP = bl.Mper.transpose() * At;   // (pairs) x NG
  MatrixXd RX = bl.MX0.transpose() * At;
  for (int i = 0; i < NG; ++i) {   // null group summed directly where its expansion does not hold
    if (gw.bad[i] || gw.W[i].null_ok) continue;
    const Weights& W = gw.W[i]; const double h = HS[i];
    for (int a = 0; a < bl.N; ++a) for (int b = a; b < bl.N; ++b) {
      const int t = pidx(a, b, bl.N);
      if (bl.kind[t] == 0) RX(bl.slot[t], i) = sp.wsum(W, 0, h, bl.col(a, 0), bl.col(b, 0));
      else if (bl.kind[t] == 1) RY(bl.slot[t], i) = sp.wsum(W, 0, h, bl.col(a, 0), bl.col(b, 0));
      else for (int s = 0; s < S; ++s) RP((size_t)s * nper + bl.slot[t], i) = sp.wsum(W, 0, h, bl.col(a, s), bl.col(b, s));
    }
  }
  auto M = [&](int i, int s, int a, int b) {
    const int t = pidx(a, b, bl.N);
    return bl.kind[t] == 0 ? RX(bl.slot[t], i) : bl.kind[t] == 1 ? RY(bl.slot[t], i) : RP((size_t)s * nper + bl.slot[t], i);
  };
  (void)ny;
  std::vector<MatrixXd> V(models.size(), MatrixXd(S, NG));
  std::vector<double> r0(p0), z(p0), u(k), c(k), a(p0 * k), T(p0 * k), Sc(k * k);
  for (int i = 0; i < NG; ++i) {
    if (gw.bad[i]) { for (auto& v : V) v.col(i).setConstant(NEG_INF); continue; }
    const double* Ai = xg.A0i[i].data(); const double ldi = gw.logdet[i];
    for (int s = 0; s < S; ++s) {
      const double yy = M(i, s, iy, iy); double q = 0;
      for (int e = 0; e < p0; ++e) r0[e] = M(i, s, e, iy);
      for (int e = 0; e < p0; ++e) { double x = 0; for (int f = 0; f < p0; ++f) x += Ai[e + f * p0] * r0[f]; z[e] = x; q += r0[e] * x; }
      for (size_t mi = 0; mi < models.size(); ++mi) {
        const std::vector<int>& mb = models[mi]; const int kk = mb.size();
        double rss = yy - q; bool full = false;
        for (int j = 0; j < kk; ++j) {
          const int cj = p0 + mb[j]; double* aj = &a[j * p0]; double* Tj = &T[j * p0];
          for (int e = 0; e < p0; ++e) aj[e] = M(i, s, e, cj);
          for (int e = 0; e < p0; ++e) { double x = 0; for (int f = 0; f < p0; ++f) x += Ai[e + f * p0] * aj[f]; Tj[e] = x; }
          c[j] = M(i, s, cj, iy);
          for (int l = 0; l <= j; ++l) Sc[j * k + l] = M(i, s, cj, p0 + mb[l]);
        }
        for (int j = 0; j < kk; ++j) {
          if (Sc[j * k + j] == 0) { for (int l = 0; l < j; ++l) Sc[j * k + l] = 0; Sc[j * k + j] = 1; u[j] = 0; continue; }   // absent column
          const double* aj = &a[j * p0]; double x = c[j];
          for (int e = 0; e < p0; ++e) x -= aj[e] * z[e];
          u[j] = x;
          for (int l = 0; l <= j; ++l) {
            if (Sc[l * k + l] == 0) continue;
            const double* Tl = &T[l * p0]; double y = 0; for (int e = 0; e < p0; ++e) y += aj[e] * Tl[e]; Sc[j * k + l] -= y;
          }
        }
        for (int j = 0; j < kk; ++j) {   // Cholesky of the Schur complement, forward solve, rss -= |L^-1 u|^2
          double dj = Sc[j * k + j]; const double dj0 = M(i, s, p0 + mb[j], p0 + mb[j]);
          for (int l = 0; l < j; ++l) dj -= Sc[j * k + l] * Sc[j * k + l];
          if (!(dj > SCHUR_PIVOT * std::max(dj0, 1e-300))) { full = true; break; }
          const double ljj = std::sqrt(dj); Sc[j * k + j] = ljj;
          for (int r = j + 1; r < kk; ++r) { double x = Sc[r * k + j]; for (int l = 0; l < j; ++l) x -= Sc[r * k + l] * Sc[j * k + l]; Sc[r * k + j] = x / ljj; }
          double x = u[j]; for (int l = 0; l < j; ++l) x -= Sc[j * k + l] * u[l]; u[j] = x / ljj; rss -= u[j] * u[j];
        }
        if (full) {   // p0's LDLT on the full system (pseudo-solve for collinear columns)
          const int p = p0 + kk; MatrixXd A(p, p); VectorXd rhs(p);
          auto cix = [&](int e) { return e < p0 ? e : p0 + mb[e - p0]; };
          for (int e = 0; e < p; ++e) { rhs[e] = M(i, s, cix(e), iy); for (int f = 0; f < p; ++f) A(e, f) = M(i, s, cix(e), cix(f)); }
          for (int j = 0; j < kk; ++j) if (A(p0 + j, p0 + j) == 0) { A(p0 + j, p0 + j) = 1; rhs[p0 + j] = 0; }
          rss = yy - rhs.dot(A.ldlt().solve(rhs));
        }
        if (!(rss > CANCEL_GRID * yy)) {   // near-perfect fit: direct residuals
          std::vector<int> cols; for (int e = 0; e < p0; ++e) cols.push_back(e);
          for (int j = 0; j < kk; ++j) cols.push_back(p0 + mb[j]);
          rss = grid_direct(bl, gw, i, s, cols);
        }
        const double sv = rss / n;
        V[mi](s, i) = sv > 0 ? -0.5 * n - 0.5 * (n * std::log(sv) + ldi) : POS_INF;
      }
    }
  }
  return V;
}

// ---------------------------------------------------------------- rotation
// U'M with repeated rows of M (samples of one individual share a genotype; the interaction column is constant
// over the reference condition): M = Z Mu, U'M = (Z'U)' Mu; only the summation order changes.
static MatrixXd rotate(const Eigen::Map<MatrixXd>& U, const MatrixXd& M) {
  const int N = M.rows(), S = M.cols();
  if (S == 0) return MatrixXd(U.cols(), 0);
  const MatrixXd Mt = M.transpose();
  std::vector<int> ord(N); for (int i = 0; i < N; ++i) ord[i] = i;
  auto less = [&](int a, int b) { const double *x = Mt.col(a).data(), *y = Mt.col(b).data();
    for (int j = 0; j < S; ++j) if (x[j] != y[j]) return x[j] < y[j]; return false; };
  std::sort(ord.begin(), ord.end(), less);
  std::vector<int> cls(N); int nc = 0;
  for (int i = 0; i < N; ++i) cls[ord[i]] = (i > 0 && !less(ord[i - 1], ord[i])) ? nc - 1 : nc++;
  if (nc > 0.75 * N) return U.transpose() * M;   // few repeats: plain product
  MatrixXd UZ = MatrixXd::Zero(U.cols(), nc), Mu(nc, S);
  for (int i = 0; i < N; ++i) { UZ.col(cls[i]) += U.row(i).transpose(); Mu.row(cls[i]) = M.row(i); }
  return UZ * Mu;
}

static bool standardize(Eigen::Ref<VectorXd> x, double& sd) {
  const int n = x.size(); x.array() -= x.mean(); sd = std::sqrt(x.squaredNorm() / (n - 1));
  if (!(sd > 0)) return false;
  x /= sd; return true;
}

// beta and se as p0: (X'W'X)^-1 by LU with W' = W / s at the optimum (A = X'WX of the fit in P)
static void beta_se(const Profile& P, int j, double& beta, double& se) {
  MatrixXd Vb = (P.A / P.s).inverse(); VectorXd bh = Vb * (P.rhs / P.s); beta = bh[j]; se = std::sqrt(Vb(j, j));
}

// ---------------------------------------------------------------- fit_gene (suez_fit_gene)
static Rcpp::List fit_gene(const Eigen::Map<MatrixXd>& U, const VectorXd& lam_, const VectorXd& Uy_, const MatrixXd& UX0_,
                           const MatrixXd& G, const MatrixXd& C, bool boot) {
  const int S = G.cols(), q = C.cols(), p0 = UX0_.cols(), n = Uy_.size();
  const ArrayXd lam = lam_.array();
  const Spectrum sp(lam); const GridW gw(sp, lam); Profile PR(sp);
  const MatrixXd UX0 = sp.permute(UX0_), Y1 = sp.permute(MatrixXd(Uy_));
  const MatrixXd MX0 = x0_moments(sp, UX0); const X0Grid xg = x0_grid(sp, gw, UX0, MX0);
  std::vector<int> c0(p0); for (int e = 0; e < p0; ++e) c0[e] = e;
  const std::vector<MatrixXd> none;

  Fit f0;
  { Block b0(sp, UX0, MX0, none, Y1, 1); MatrixXd V0 = grid_values(b0, gw, xg, std::vector<std::vector<int>>(1))[0];
    f0 = refine(PR, make_model(b0, 0, c0), V0.data()); }
  if (!f0.ok) Rcpp::stop("no-genotype model did not converge");

  MatrixXd Gs = G; VectorXd gsd(S);
  for (int s = 0; s < S; ++s) standardize(Gs.col(s), gsd[s]);
  std::vector<MatrixXd> Is(q); std::vector<std::vector<bool>> iok(q, std::vector<bool>(S)); VectorXi df = VectorXi::Zero(S);
  for (int j = 0; j < q; ++j) {
    Is[j] = G.array().colwise() * C.col(j).array();
    for (int s = 0; s < S; ++s) { double sd; iok[j][s] = standardize(Is[j].col(s), sd); if (!iok[j][s]) Is[j].col(s).setZero(); else ++df[s]; }
  }
  std::vector<int> all(1 + q); for (int j = 0; j <= q; ++j) all[j] = j;
  const std::vector<std::vector<int>> mdata = q ? std::vector<std::vector<int>>{{0}, all} : std::vector<std::vector<int>>{{0}};
  const std::vector<std::vector<int>> mboot = q ? std::vector<std::vector<int>>{{}, {0}, all} : std::vector<std::vector<int>>{{}, {0}};
  auto cols1 = [&]() { std::vector<int> c = c0; c.push_back(p0); return c; };
  auto cols2 = [&](int s) { std::vector<int> c = cols1(); for (int j = 0; j < q; ++j) if (iok[j][s]) c.push_back(p0 + 1 + j); return c; };
  auto ok_nested = [](Fit& f, const Fit& g0) { if (!f.ok || f.value < g0.value - 1e-6) f.ok = false; };

  Rcpp::NumericMatrix out(S, 7);   // l_geno, l_interact, df, beta, se, l_boot_geno, l_boot_interact
  std::fill(out.begin(), out.end(), NA_REAL);
  Rcpp::RNGScope rng;
  std::vector<int> pos(n); for (int i = 0; i < n; ++i) pos[sp.perm[i]] = i;
  for (int s0 = 0; s0 < S; s0 += SNP_BLOCK) {
    const int sb = std::min(SNP_BLOCK, S - s0);
    std::vector<MatrixXd> B(1, sp.permute(rotate(U, Gs.middleCols(s0, sb))));
    for (int j = 0; j < q; ++j) B.push_back(sp.permute(rotate(U, Is[j].middleCols(s0, sb))));
    std::vector<Fit> f1(sb);
    {
      Block bd(sp, UX0, MX0, B, Y1, sb);
      std::vector<MatrixXd> VV = grid_values(bd, gw, xg, mdata);
      for (int sl = 0; sl < sb; ++sl) {
        const int s = s0 + sl; out(s, 2) = df[s];
        VectorXd v1 = VV[0].row(sl).transpose(); f1[sl] = refine(PR, make_model(bd, sl, cols1()), v1.data()); ok_nested(f1[sl], f0);
        if (!f1[sl].ok) continue;
        out(s, 0) = f1[sl].value;
        double be, se; beta_se(PR, p0, be, se); out(s, 3) = be / gsd[s]; out(s, 4) = se / gsd[s];
        if (df[s] > 0) {
          VectorXd v2 = VV[1].row(sl).transpose(); Fit f2 = refine(PR, make_model(bd, sl, cols2(s)), v2.data()); ok_nested(f2, f1[sl]);
          if (f2.ok) out(s, 1) = f2.value;
        }
      }
    }
    if (!boot) continue;
    // bootstrap: one y per SNP, drawn in SNP order and in the original contrast order (as p0)
    MatrixXd Yb = MatrixXd::Zero(n, sb);
    for (int sl = 0; sl < sb; ++sl) {
      const Fit& f = f1[sl]; if (!f.ok) continue;
      for (int i = 0; i < n; ++i) Yb(pos[i], sl) = std::sqrt(f.s2k * lam[i] + f.s2) * R::norm_rand();
      Yb.col(sl) += UX0 * f.b.head(p0) + B[0].col(sl) * f.b[p0];
    }
    Block bb(sp, UX0, MX0, B, Yb, sb);
    std::vector<MatrixXd> WW = grid_values(bb, gw, xg, mboot);
    for (int sl = 0; sl < sb; ++sl) {
      const int s = s0 + sl; if (!f1[sl].ok) continue;
      VectorXd v = WW[0].row(sl).transpose(); Fit g0 = refine(PR, make_model(bb, sl, c0), v.data());
      if (!g0.ok) continue;
      v = WW[1].row(sl).transpose(); Fit g1 = refine(PR, make_model(bb, sl, cols1()), v.data()); ok_nested(g1, g0);
      if (!g1.ok) continue;
      out(s, 5) = g1.value;
      if (df[s] > 0) {
        v = WW[2].row(sl).transpose(); Fit g2 = refine(PR, make_model(bb, sl, cols2(s)), v.data()); ok_nested(g2, g1);
        if (g2.ok) out(s, 6) = g2.value;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("l0") = f0.value, Rcpp::Named("res") = out);
}

// [[Rcpp::export]]
Rcpp::List fit_gene_proto(Eigen::Map<Eigen::MatrixXd> U, Eigen::VectorXd lambda, Eigen::VectorXd Uy, Eigen::MatrixXd UX0,
                          Eigen::MatrixXd G, Eigen::MatrixXd C, bool boot) { return fit_gene(U, lambda, Uy, UX0, G, C, boot); }
