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
static const int MAX_PEAKS = 3;
static long FALLBACKS = 0;
// [[Rcpp::export]]
double fallbacks() { double f = FALLBACKS; FALLBACKS = 0; return f; }   // local maxima of the grid refined per fit

struct Fit { double value = NAN, s2k = NAN, s2 = NAN; VectorXd b; bool ok = false; };

// Profile log-likelihood at h (no 2 pi constant) and its first two derivatives in h:
//   f   = -(n/2) log(g) - (1/2) sum log d,      d = h lam + 1 - h,  c = lam - 1,  w = 1/d,  g = y'Py
//   f'  = (n/2) Q1/g - (1/2) sum c w,            Q1 = sum c (Py)^2
//   f'' = -(n/2) (2 Q2/g - Q1^2/g^2) + (1/2) sum c^2 w^2,   Q2 = v'Pv, v = c * Py
// with P = W - W X (X'WX)^-1 X'W, so Py = W r for the GLS residual r. Sets b and s (= g / n).
// p1: same formulas, no per-call heap allocation (buffers reused across calls), no Xw copy, c = lam - 1
// computed once per gene.
struct Profile {
  ArrayXd lam, c, c2, d, w, r, Py, v; MatrixXd A; VectorXd b, rhs, Xv, tmp; Eigen::LDLT<MatrixXd> L; double s = NAN;
  explicit Profile(const ArrayXd& lam_) : lam(lam_), c(lam_ - 1), c2(c.square()) {}
  double operator()(const MatrixXd& X, const VectorXd& y, const ArrayXd&, double h, double* f1 = nullptr, double* f2 = nullptr) {
    const int n = y.size(), p = X.cols();
    d = h * lam + (1 - h);
    if ((d <= 0).any()) { if (f1) *f1 = *f2 = NAN; return NEG_INF; }
    w = d.inverse(); r = y.array();
    if (p) {
      A.resize(p, p); rhs.resize(p);
      for (int a = 0; a < p; ++a) {
        Py = X.col(a).array() * w;                     // scratch: w x_a
        rhs[a] = (Py * y.array()).sum();
        for (int e = 0; e <= a; ++e) A(a, e) = A(e, a) = (Py * X.col(e).array()).sum();
      }
      L.compute(A); b = L.solve(rhs); r.matrix().noalias() -= X * b;
    } else b.resize(0);
    Py = w * r; double g = (r * Py).sum(); s = g / n;
    if (!(s > 0)) { if (f1) *f1 = *f2 = NAN; return std::numeric_limits<double>::infinity(); }
    if (f1) {
      v = c * Py;
      double Q1 = (v * Py).sum(), Q2 = (w * v.square()).sum();
      if (p) { Xv.noalias() = X.transpose() * (w * v).matrix(); tmp = L.solve(Xv); Q2 -= Xv.dot(tmp); }
      *f1 = 0.5 * n * Q1 / g - 0.5 * (c * w).sum();
      *f2 = -0.5 * n * (2 * Q2 / g - Q1 * Q1 / (g * g)) + 0.5 * (c2 * w.square()).sum();
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
static Fit refine(Profile& P, const MatrixXd& X, const VectorXd& y, const ArrayXd& lam, const double* v, int best) {
  std::vector<int> peaks;
  for (int i = 0; i < NG; ++i)
    if (std::isfinite(v[i]) && (i == 0 || v[i] >= v[i - 1]) && (i == NG - 1 || v[i] >= v[i + 1])) peaks.push_back(i);
  std::sort(peaks.begin(), peaks.end(), [&](int p, int q) { return v[p] > v[q]; });
  if (peaks.empty()) peaks.push_back(best);
  double hb = HS[peaks[0]], vb = NEG_INF;
  for (size_t j = 0; j < peaks.size() && j < (size_t)MAX_PEAKS; ++j) {
    auto m = newton(P, X, y, lam, v, peaks[j]);
    if (m.second > vb) { vb = m.second; hb = m.first; }
  }
  Fit f; f.value = P(X, y, lam, hb); f.s2k = P.s * hb; f.s2 = P.s * (1 - hb); f.b = P.b;
  f.ok = std::isfinite(f.value);
  return f;
}

// p2: grid values as GEMMs. Every quantity grid_values needs is sum_t w_t(h_i) u_t v_t for a pair of columns
// (u, v) of [X0, B_1..B_k, Y]: stacking the products u .* v as columns of P (n x m), all grid points at once
// are Wg * P with Wg(i, t) = 1 / (h_i lam_t + 1 - h_i) (NG x n). Per SNP and grid point the normal
// equations are then solved by the Schur complement on the shared X0 block (A00^-1 once per grid point):
//   rss = yy - r0' A00^-1 r0 - u' Sc^-1 u,  T = A00^-1 a,  Sc = Bb - a'T,  u = c - T'r0.
// The nested models (first kk blocks, for each kk in kks) share one set of moments. SNPs go in chunks of SB so
// the product buffer stays small. A (near-)singular Schur pivot falls back to the baseline LDLT on full A.
struct Grid {
  MatrixXd W; VectorXd logdet; std::vector<char> bad;
  explicit Grid(const ArrayXd& lam) : W(NG, lam.size()), logdet(NG), bad(NG) {
    for (int i = 0; i < NG; ++i) {
      ArrayXd d = HS[i] * lam + (1 - HS[i]);
      bad[i] = (d <= 0).any();
      if (bad[i]) { W.row(i).setZero(); logdet[i] = NAN; continue; }
      W.row(i) = d.inverse().matrix().transpose(); logdet[i] = d.log().sum();
    }
  }
};
static const int SB = 64;

static std::vector<MatrixXd> grid_values_multi(const Grid& gr, const MatrixXd& X0, const std::vector<MatrixXd>& B,
                                               const MatrixXd& Y, int S, const std::vector<int>& kks) {
  const int n = X0.rows(), p0 = X0.cols(), k = B.size(), N = p0 + k + 1, iy = p0 + k;
  const bool shared_y = Y.cols() == 1;
  auto is_shared = [&](int a) { return a < p0 || (a == iy && shared_y); };
  // term index of each pair (a <= b): shared ones into Psh, per-SNP ones into Pc
  std::vector<int> idx(N * N, -1); std::vector<std::pair<int, int>> sh, ps;
  for (int a = 0; a < N; ++a) for (int b = a; b < N; ++b) {
    if (is_shared(a) && is_shared(b)) { idx[a * N + b] = idx[b * N + a] = sh.size(); sh.push_back({a, b}); }
    else { idx[a * N + b] = idx[b * N + a] = ps.size(); ps.push_back({a, b}); }
  }
  auto col = [&](int a, int s) -> Eigen::Ref<const VectorXd> {
    return a < p0 ? X0.col(a) : a < iy ? B[a - p0].col(s) : Y.col(shared_y ? 0 : s);
  };
  MatrixXd Psh(n, sh.size());
  for (size_t t = 0; t < sh.size(); ++t) Psh.col(t) = col(sh[t].first, 0).cwiseProduct(col(sh[t].second, 0));
  MatrixXd Rsh = gr.W * Psh;                                       // NG x shared terms
  // per grid point: A00^-1, and with a shared y, r0 and q0 = r0' A00^-1 r0
  std::vector<MatrixXd> A0i(NG); MatrixXd z0(p0, NG); VectorXd q0(NG);
  for (int i = 0; i < NG; ++i) {
    if (gr.bad[i]) continue;
    MatrixXd A00(p0, p0);
    for (int a = 0; a < p0; ++a) for (int b = 0; b < p0; ++b) A00(a, b) = Rsh(i, idx[a * N + b]);
    A0i[i] = A00.ldlt().solve(MatrixXd::Identity(p0, p0));
    if (shared_y) {
      VectorXd r0(p0); for (int a = 0; a < p0; ++a) r0[a] = Rsh(i, idx[a * N + iy]);
      z0.col(i) = A0i[i] * r0; q0[i] = r0.dot(z0.col(i));
    }
  }
  std::vector<MatrixXd> V(kks.size(), MatrixXd(S, NG));
  const int m = ps.size();
  MatrixXd P(n, m * SB), R(NG, m * SB);
  VectorXd r0(p0), z(p0), u(k), c(k); MatrixXd a(p0, k), T(p0, k), Sc(k, k);
  for (int s0 = 0; s0 < S; s0 += SB) {
    const int sb = std::min(SB, S - s0);
    for (int t = 0; t < m; ++t) for (int sl = 0; sl < sb; ++sl)
      P.col(t * sb + sl) = col(ps[t].first, s0 + sl).cwiseProduct(col(ps[t].second, s0 + sl));
    R.leftCols(m * sb).noalias() = gr.W * P.leftCols(m * sb);
    // moment of pair (a, b) for SNP sl at grid point i
    auto M = [&](int i, int sl, int a, int b) { int t = idx[a * N + b]; return (is_shared(a) && is_shared(b)) ? Rsh(i, t) : R(i, t * sb + sl); };
    for (int i = 0; i < NG; ++i) {
      if (gr.bad[i]) { for (auto& v : V) v.block(s0, i, sb, 1).setConstant(NEG_INF); continue; }
      for (int sl = 0; sl < sb; ++sl) {
        double yy = M(i, sl, iy, iy), q = 0;
        if (shared_y) { q = q0[i]; z = z0.col(i); }
        else if (p0) { for (int e = 0; e < p0; ++e) r0[e] = M(i, sl, e, iy); z.noalias() = A0i[i] * r0; q = r0.dot(z); }
        for (size_t mi = 0; mi < kks.size(); ++mi) {
          const int kk = kks[mi]; double rss = yy - q; bool fb = false;
          if (kk) {
            for (int j = 0; j < kk; ++j) {
              for (int e = 0; e < p0; ++e) a(e, j) = M(i, sl, e, p0 + j);
              c[j] = M(i, sl, p0 + j, iy);
              for (int l = 0; l <= j; ++l) Sc(j, l) = M(i, sl, p0 + j, p0 + l);
            }
            auto aa = a.leftCols(kk); auto TT = T.leftCols(kk);
            if (p0) TT.noalias() = A0i[i] * aa;
            for (int j = 0; j < kk; ++j) {
              if (Sc(j, j) == 0) { for (int l = 0; l < j; ++l) Sc(j, l) = 0; Sc(j, j) = 1; u[j] = 0; continue; }   // absent column
              u[j] = c[j] - (p0 ? aa.col(j).dot(z) : 0.0);
              for (int l = 0; l <= j; ++l) if (Sc(l, l) != 0 || l == j) Sc(j, l) -= p0 ? aa.col(j).dot(TT.col(l)) : 0.0;
            }
            // in-place Cholesky of the kk x kk Schur complement (lower), then rss -= |L^-1 u|^2
            for (int j = 0; j < kk && !fb; ++j) {
              double dj = Sc(j, j); const double dj0 = dj;
              for (int l = 0; l < j; ++l) dj -= Sc(j, l) * Sc(j, l);
              if (!(dj > 1e-10 * std::max(dj0, 1e-300))) { fb = true; break; }
              Sc(j, j) = std::sqrt(dj);
              for (int r = j + 1; r < kk; ++r) { double x = Sc(r, j); for (int l = 0; l < j; ++l) x -= Sc(r, l) * Sc(j, l); Sc(r, j) = x / Sc(j, j); }
              double x = u[j]; for (int l = 0; l < j; ++l) x -= Sc(j, l) * u[l]; u[j] = x / Sc(j, j); rss -= u[j] * u[j];
            }
            if (fb) {   // baseline: LDLT on the full (p0 + kk) system
              const int p = p0 + kk; MatrixXd A(p, p); VectorXd rhs(p);
              for (int e = 0; e < p; ++e) { rhs[e] = M(i, sl, e, iy); for (int f = 0; f < p; ++f) A(e, f) = M(i, sl, e, f); }
              for (int j = 0; j < kk; ++j) if (A(p0 + j, p0 + j) == 0) { A(p0 + j, p0 + j) = 1; rhs[p0 + j] = 0; }
              rss = yy - rhs.dot(A.ldlt().solve(rhs)); ++FALLBACKS;
            }
          }
          double sv = rss / n;
          V[mi](s0 + sl, i) = sv > 0 ? -0.5 * n - 0.5 * (n * std::log(sv) + gr.logdet[i]) : std::numeric_limits<double>::infinity();
        }
      }
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
static Rcpp::List fit_gene(const Eigen::Map<MatrixXd>& U, const VectorXd& lam_, const VectorXd& Uy, const MatrixXd& UX0,
                           const MatrixXd& G, const MatrixXd& C, bool boot) {
  const int S = G.cols(), q = C.cols(), p0 = UX0.cols(), n = Uy.size();
  const ArrayXd lam = lam_.array();
  MatrixXd Y1(n, 1); Y1.col(0) = Uy;
  Profile PR(lam);

  Grid gr(lam);
  MatrixXd V0 = grid_values_multi(gr, UX0, {}, Y1, 1, {0})[0];
  Fit f0 = refine(PR, UX0, Uy, lam, V0.data(), argmax(V0.data()));
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
  std::vector<MatrixXd> VV = grid_values_multi(gr, UX0, q ? blocks : b1, Y1, S, q ? std::vector<int>{1, 1 + q} : std::vector<int>{1});
  const MatrixXd& V1 = VV[0]; MatrixXd V2 = q ? VV[1] : MatrixXd();
  std::vector<Fit> f1(S);
  for (int s = 0; s < S; ++s) {
    out(s, 2) = df[s];
    MatrixXd X1 = X1_of(s);
    VectorXd v1 = V1.row(s).transpose(); f1[s] = refine(PR, X1, Uy, lam, v1.data(), argmax(v1.data())); ok_nested(f1[s], f0);
    if (!f1[s].ok) continue;
    out(s, 0) = f1[s].value;
    out(s, 3) = f1[s].b[p0] / gsd[s]; out(s, 4) = std::sqrt(PR.s * PR.L.solve(VectorXd::Unit(p0 + 1, p0))[p0]) / gsd[s];   // PR holds the fit at the optimum: Vb = s A^-1
    if (df[s] > 0) {
      VectorXd v2 = V2.row(s).transpose(); Fit f2 = refine(PR, X2_of(s), Uy, lam, v2.data(), argmax(v2.data())); ok_nested(f2, f1[s]);
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
  std::vector<MatrixXd> WW = grid_values_multi(gr, UX0, q ? blocks : b1, Yb, S, q ? std::vector<int>{0, 1, 1 + q} : std::vector<int>{0, 1});
  const MatrixXd &W0 = WW[0], &W1 = WW[1]; MatrixXd W2 = q ? WW[2] : MatrixXd();
  for (int s = 0; s < S; ++s) {
    if (!f1[s].ok) continue;
    VectorXd y = Yb.col(s), v;
    v = W0.row(s).transpose(); Fit g0 = refine(PR, UX0, y, lam, v.data(), argmax(v.data()));
    if (!g0.ok) continue;
    v = W1.row(s).transpose(); Fit g1 = refine(PR, X1_of(s), y, lam, v.data(), argmax(v.data())); ok_nested(g1, g0);
    if (!g1.ok) continue;
    out(s, 5) = g1.value;
    if (df[s] > 0) {
      v = W2.row(s).transpose(); Fit g2 = refine(PR, X2_of(s), y, lam, v.data(), argmax(v.data())); ok_nested(g2, g1);
      if (g2.ok) out(s, 6) = g2.value;
    }
  }
  return Rcpp::List::create(Rcpp::Named("l0") = f0.value, Rcpp::Named("res") = out);
}


// [[Rcpp::export]]
Rcpp::List fit_gene_proto(Eigen::Map<Eigen::MatrixXd> U, Eigen::VectorXd lambda, Eigen::VectorXd Uy, Eigen::MatrixXd UX0,
                          Eigen::MatrixXd G, Eigen::MatrixXd C, bool boot) { return fit_gene(U, lambda, Uy, UX0, G, C, boot); }
