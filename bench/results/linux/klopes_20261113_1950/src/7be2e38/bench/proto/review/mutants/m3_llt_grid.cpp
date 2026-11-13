// [[Rcpp::depends(RcppEigen)]]
// Exact ML fits for map_interaction_qtl, all cis-SNPs of one gene per call.
// Model on rotated data: Uy = UX b + e, e ~ N(0, diag(sigma2_k * lambda + sigma2)). Writing
// sigma2_k = s h, sigma2 = s (1 - h), b and s have closed forms for each h in [0, 1], so each fit is a 1-D
// maximisation over h. The search follows GEMMA / FaST-LMM: a grid uniform in x = log10(delta),
// delta = sigma2_k / sigma2 = h / (1 - h), evaluated for every SNP of the gene together (shared condition
// terms once, genotype terms as matrix products across SNPs); then safeguarded Newton steps with analytic
// derivatives around every local maximum of the grid (the likelihood can be bimodal), keeping the best.
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <limits>
using Eigen::ArrayXd; using Eigen::ArrayXXd; using Eigen::MatrixXd; using Eigen::VectorXd; using Eigen::VectorXi;
static const double NEG_INF = -std::numeric_limits<double>::infinity();
static const double L10 = 2.302585092994046;   // log(10)

// h grid: h = 0, then x = log10(delta) uniform on [-5, 5] (80 points, spacing 0.13), then h = 1. Uniform in
// x, the grid suits kernels of any scale; in h a fixed grid misses narrow peaks when the kernel's
// eigenvalues span orders of magnitude (e.g. a learned cov_suez). Peaks beyond |x| = 5 are reached from
// the end points (Newton works up to |x| = 12).
static std::vector<double> make_grid() {
  const int nx = 80; const double xl = 5;
  std::vector<double> h(1, 0.0);
  for (int i = 0; i < nx; ++i) h.push_back(1 / (1 + std::pow(10.0, -(-xl + 2 * xl * i / (nx - 1.0)))));
  h.push_back(1); return h;
}
static const std::vector<double> HS = make_grid();
static const int NG = HS.size();
static const int MAX_PEAKS = 3;   // local maxima of the grid refined per fit

struct Fit { double value = NAN, s2k = NAN, s2 = NAN; VectorXd b; bool ok = false; };

// Profile log-likelihood at h (no 2 pi constant) and its first two derivatives in h:
//   f   = -(n/2) log(g) - (1/2) sum log d,      d = h lam + 1 - h,  c = lam - 1,  w = 1/d,  g = y'Py
//   f'  = (n/2) Q1/g - (1/2) sum c w,            Q1 = sum c (Py)^2
//   f'' = -(n/2) (2 Q2/g - Q1^2/g^2) + (1/2) sum c^2 w^2,   Q2 = v'Pv, v = c * Py
// with P = W - W X (X'WX)^-1 X'W, so Py = W r for the GLS residual r. Sets b and s (= g / n).
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

// Safeguarded Newton in x = log10(delta) on the bracket around grid point `best`: h = 1 / (1 + 10^-x),
// dh/dx = ln10 h (1 - h), d2h/dx2 = ln10 (1 - 2h) dh/dx. Steps leaving the bracket become bisections.
// Returns the best h found and its value.
static std::pair<double, double> newton(Profile& P, const MatrixXd& X, const VectorXd& y, const ArrayXd& lam,
                                        const double* v, int best) {
  const double XCAP = 12;
  auto X_of = [&](double h) { return h <= 0 ? -XCAP : h >= 1 ? XCAP : std::max(-XCAP, std::min(XCAP, std::log10(h / (1 - h)))); };
  auto H = [](double x) { return 1 / (1 + std::pow(10.0, -x)); };
  double hb = HS[best], vb = v[best], f1, f2;
  double a = HS[std::max(0, best - 1)], c = HS[std::min(NG - 1, best + 1)];
  if (!std::isfinite(vb) || !(c > a)) return {hb, vb};
  double lo = X_of(a), hi = X_of(c), x = X_of(hb);
  for (int it = 0; it < 60; ++it) {
    double h = H(x), fv = P(X, y, lam, h, &f1, &f2);
    if (fv > vb) { vb = fv; hb = h; }
    double dh = L10 * h * (1 - h), d2h = L10 * (1 - 2 * h) * dh;
    double g1 = f1 * dh, g2 = f2 * dh * dh + f1 * d2h;
    if (!std::isfinite(g1)) break;
    if (g1 > 0) lo = std::max(lo, x); else hi = std::min(hi, x);       // the maximum is uphill
    double xn = g2 < 0 ? x - g1 / g2 : NAN;
    if (!std::isfinite(xn) || xn <= lo || xn >= hi) xn = 0.5 * (lo + hi);
    if (std::fabs(xn - x) < 1e-9 || hi - lo < 1e-10) {
      double hh = H(xn), fv2 = P(X, y, lam, hh);
      if (fv2 > vb) { vb = fv2; hb = hh; }
      break;
    }
    x = xn;
  }
  return {hb, vb};
}

// Refine every local maximum of the grid values v (up to MAX_PEAKS, highest first) and return the best fit.
static Fit refine(const MatrixXd& X, const VectorXd& y, const ArrayXd& lam, const double* v, int best) {
  std::vector<int> peaks;
  for (int i = 0; i < NG; ++i)
    if (std::isfinite(v[i]) && (i == 0 || v[i] >= v[i - 1]) && (i == NG - 1 || v[i] >= v[i + 1])) peaks.push_back(i);
  std::sort(peaks.begin(), peaks.end(), [&](int p, int q) { return v[p] > v[q]; });
  if (peaks.empty()) peaks.push_back(best);
  Profile P; double hb = HS[peaks[0]], vb = NEG_INF;
  for (size_t j = 0; j < peaks.size() && j < (size_t)MAX_PEAKS; ++j) {
    auto m = newton(P, X, y, lam, v, peaks[j]);
    if (m.second > vb) { vb = m.second; hb = m.first; }
  }
  Fit f; f.value = P(X, y, lam, hb); f.s2k = P.s * hb; f.s2 = P.s * (1 - hb); f.b = P.b;
  f.ok = std::isfinite(f.value);
  return f;
}

// Grid values V (S x NG) for models X_s = [X0, B_1[,s], ..., B_k[,s]], y_s = Y[,s] (or Y[,0] if one column).
// A zero column in a block (constant interaction) is left out of that SNP's model.
static MatrixXd grid_values(const MatrixXd& X0, const std::vector<MatrixXd>& B, const MatrixXd& Y, const ArrayXd& lam, int S) {
  const int n = X0.rows(), p0 = X0.cols(), k = B.size(), p = p0 + k;
  const bool shared_y = Y.cols() == 1;
  MatrixXd V(S, NG);
  std::vector<VectorXd> bb(k * k);                       // colsum(B_j .* W .* B_l)
  std::vector<MatrixXd> A0B(k);                           // X0' W B_j   (p0 x S)
  std::vector<VectorXd> By(k);                            // B_j' W y_s  (S)
  for (int i = 0; i < NG; ++i) {
    double h = HS[i];
    ArrayXd d = h * lam + (1 - h);
    if ((d <= 0).any()) { V.col(i).setConstant(NEG_INF); continue; }
    ArrayXd w = d.inverse(); double logdet = d.log().sum();
    MatrixXd X0w = (X0.array().colwise() * w).matrix();
    MatrixXd A00 = X0w.transpose() * X0;
    MatrixXd r0 = X0w.transpose() * Y;                   // p0 x (1 or S)
    ArrayXd yy = (Y.array().square().colwise() * w).colwise().sum().transpose();
    std::vector<ArrayXXd> Bw(k);
    for (int j = 0; j < k; ++j) {
      Bw[j] = B[j].array().colwise() * w;
      A0B[j].noalias() = X0w.transpose() * B[j];
      By[j] = shared_y ? VectorXd(Bw[j].matrix().transpose() * Y.col(0))
                       : VectorXd((Bw[j] * Y.array()).colwise().sum().transpose());
      for (int l = 0; l <= j; ++l) bb[j * k + l] = (Bw[j] * B[l].array()).colwise().sum().transpose();
    }
    MatrixXd A(p, p); VectorXd rhs(p);
    for (int s = 0; s < S; ++s) {
      A.topLeftCorner(p0, p0) = A00;
      rhs.head(p0) = r0.col(shared_y ? 0 : s);
      for (int j = 0; j < k; ++j) {
        A.block(0, p0 + j, p0, 1) = A0B[j].col(s); A.block(p0 + j, 0, 1, p0) = A0B[j].col(s).transpose();
        for (int l = 0; l <= j; ++l) A(p0 + j, p0 + l) = A(p0 + l, p0 + j) = bb[j * k + l][s];
        rhs[p0 + j] = By[j][s];
        if (bb[j * k + j][s] == 0) { A(p0 + j, p0 + j) = 1; rhs[p0 + j] = 0; }   // absent column
      }
      VectorXd b = A.llt().solve(rhs);
      double rss = yy[shared_y ? 0 : s] - rhs.dot(b), sv = rss / n;
      V(s, i) = sv > 0 ? -0.5 * n - 0.5 * (n * std::log(sv) + logdet) : std::numeric_limits<double>::infinity();
    }
  }
  return V;
}

static int argmax(const double* v) { int b = 0; for (int i = 1; i < NG; ++i) if (v[i] > v[b]) b = i; return b; }

static bool standardize(Eigen::Ref<VectorXd> x, double& sd) {
  const int n = x.size(); x.array() -= x.mean(); sd = std::sqrt(x.squaredNorm() / (n - 1));
  if (!(sd > 0)) return false;
  x /= sd; return true;
}

// U: samples x contrasts; G: samples x SNPs (genotypes); C: samples x q interaction multipliers (the
// numeric condition, or indicators of the non-reference levels). Returns l0 and an S x 7 matrix:
// l_geno, l_interact, df_interact, beta_geno, se_geno (per allele), l_boot_geno, l_boot_interact.
static Rcpp::List fit_gene(const MatrixXd& U, const VectorXd& lam_, const VectorXd& Uy, const MatrixXd& UX0,
                           const MatrixXd& G, const MatrixXd& C, bool boot) {
  const int S = G.cols(), q = C.cols(), p0 = UX0.cols(), n = Uy.size();
  const ArrayXd lam = lam_.array();
  MatrixXd Y1(n, 1); Y1.col(0) = Uy;

  MatrixXd V0 = grid_values(UX0, {}, Y1, lam, 1);
  Fit f0 = refine(UX0, Uy, lam, V0.data(), argmax(V0.data()));
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
  const MatrixXd& UG = blocks[0];
  auto X1_of = [&](int s) { MatrixXd X(n, p0 + 1); X << UX0, UG.col(s); return X; };
  auto X2_of = [&](int s) {
    MatrixXd X(n, p0 + 1 + df[s]); X.leftCols(p0) = UX0; X.col(p0) = UG.col(s);
    for (int j = 0, c = p0 + 1; j < q; ++j) if (iok[j][s]) X.col(c++) = blocks[1 + j].col(s);
    return X;
  };
  auto ok_nested = [](Fit& f, const Fit& g0) { if (!f.ok || f.value < g0.value - 1e-6) f.ok = false; };

  Rcpp::NumericMatrix out(S, 7);   // l_geno, l_interact, df, beta, se, l_boot_geno, l_boot_interact
  std::fill(out.begin(), out.end(), NA_REAL);
  std::vector<MatrixXd> b1(blocks.begin(), blocks.begin() + 1);

  // genotype and interaction models on the data, all SNPs per grid point
  MatrixXd V1 = grid_values(UX0, b1, Y1, lam, S);
  MatrixXd V2 = q ? grid_values(UX0, blocks, Y1, lam, S) : MatrixXd();
  std::vector<Fit> f1(S);
  for (int s = 0; s < S; ++s) {
    out(s, 2) = df[s];
    MatrixXd X1 = X1_of(s);
    VectorXd v1 = V1.row(s).transpose(); f1[s] = refine(X1, Uy, lam, v1.data(), argmax(v1.data())); ok_nested(f1[s], f0);
    if (!f1[s].ok) continue;
    out(s, 0) = f1[s].value;
    ArrayXd w = (f1[s].s2k * lam + f1[s].s2).inverse();
    MatrixXd Xw = (X1.array().colwise() * w).matrix();
    MatrixXd Vb = (Xw.transpose() * X1).inverse();
    VectorXd bh = Vb * (Xw.transpose() * Uy);
    out(s, 3) = bh[p0] / gsd[s]; out(s, 4) = std::sqrt(Vb(p0, p0)) / gsd[s];
    if (df[s] > 0) {
      VectorXd v2 = V2.row(s).transpose(); Fit f2 = refine(X2_of(s), Uy, lam, v2.data(), argmax(v2.data())); ok_nested(f2, f1[s]);
      if (f2.ok) out(s, 1) = f2.value;
    }
  }
  if (!boot) return Rcpp::List::create(Rcpp::Named("l0") = f0.value, Rcpp::Named("res") = out);

  // bootstrap: one simulated y per SNP, drawn in SNP order like the R code
  Rcpp::RNGScope rng;
  MatrixXd Yb = MatrixXd::Zero(n, S);
  for (int s = 0; s < S; ++s) {
    if (!f1[s].ok) continue;
    for (int i = 0; i < n; ++i) Yb(i, s) = std::sqrt(f1[s].s2k * lam[i] + f1[s].s2) * R::norm_rand();
    Yb.col(s) += X1_of(s) * f1[s].b;
  }
  MatrixXd W0 = grid_values(UX0, {}, Yb, lam, S), W1 = grid_values(UX0, b1, Yb, lam, S);
  MatrixXd W2 = q ? grid_values(UX0, blocks, Yb, lam, S) : MatrixXd();
  for (int s = 0; s < S; ++s) {
    if (!f1[s].ok) continue;
    VectorXd y = Yb.col(s), v;
    v = W0.row(s).transpose(); Fit g0 = refine(UX0, y, lam, v.data(), argmax(v.data()));
    if (!g0.ok) continue;
    v = W1.row(s).transpose(); Fit g1 = refine(X1_of(s), y, lam, v.data(), argmax(v.data())); ok_nested(g1, g0);
    if (!g1.ok) continue;
    out(s, 5) = g1.value;
    if (df[s] > 0) {
      v = W2.row(s).transpose(); Fit g2 = refine(X2_of(s), y, lam, v.data(), argmax(v.data())); ok_nested(g2, g1);
      if (g2.ok) out(s, 6) = g2.value;
    }
  }
  return Rcpp::List::create(Rcpp::Named("l0") = f0.value, Rcpp::Named("res") = out);
}


// [[Rcpp::export]]
Rcpp::List fit_gene_proto(Eigen::MatrixXd U, Eigen::VectorXd lambda, Eigen::VectorXd Uy, Eigen::MatrixXd UX0,
                          Eigen::MatrixXd G, Eigen::MatrixXd C, bool boot) { return fit_gene(U, lambda, Uy, UX0, G, C, boot); }
